"""Readers for the 1D NMR formats this lab produces.

Supported inputs
----------------
* Bruker TopSpin experiment folders, either on disk or inside a ``.zip``
  (raw ``fid`` plus any processed ``pdata/N`` data).
* ACD/Labs ``.esp`` spectra exported from TopSpin.
* JCAMP-DX ``.jdx`` / ``.dx`` files.

Every reader returns :class:`Spectrum` objects that share one convention,
copied from Bruker: **index 0 is the left-hand (high ppm) edge**.
"""

from __future__ import annotations

import math
import os
import re
import struct
import zipfile

from . import dsp


class Spectrum:
    """One 1D spectrum plus whatever raw data came with it."""

    def __init__(self, name, real, sf, sw_hz, offset_ppm, imag=None,
                 fid=None, meta=None, source=""):
        self.name = name
        self.real = real                    # processed, high ppm first
        self.imag = imag
        self.fid = fid                      # complex time domain, or None
        self.sf = sf                        # observe frequency, MHz
        self.sw_hz = sw_hz                  # spectral width, Hz
        self.offset_ppm = offset_ppm        # ppm of real[0]
        self.meta = meta or {}
        self.source = source
        self.source_index = 0       # position within the file it was read from

        # Everything below is user state that the GUI mutates.
        self.p0 = 0.0
        self.p1 = 0.0
        self.lb = float(self.meta.get("LB", 0.0) or 0.0)
        self.si = len(real)
        self.ref_shift = 0.0                # manual ppm calibration offset
        self.baseline_on = False
        self.baseline_order = 3
        self.baseline_method = "spline"
        self.visible = True
        self.color = "#1f77b4"
        self.scale = 1.0
        self.v_offset = 0.0

        # Analysis results owned by the GUI.
        self.regions = []
        self.peaks = []

    # -- axis ---------------------------------------------------------------

    @property
    def npoints(self):
        return len(self.real)

    @property
    def sw_ppm(self):
        return self.sw_hz / self.sf if self.sf else 0.0

    @property
    def delta_ppm(self):
        """Signed ppm step between adjacent points (negative: ppm decreases)."""
        n = self.npoints
        return -self.sw_ppm / n if n else 0.0

    def ppm(self, index):
        return self.offset_ppm + self.ref_shift + index * self.delta_ppm

    def index(self, ppm):
        d = self.delta_ppm
        if d == 0.0:
            return 0
        return int(round((ppm - self.offset_ppm - self.ref_shift) / d))

    def clamp(self, index):
        return max(0, min(self.npoints - 1, index))

    @property
    def limits(self):
        """(left_ppm, right_ppm) of the full spectrum."""
        return self.ppm(0), self.ppm(self.npoints - 1)

    def hz_per_point(self):
        return self.sw_hz / self.npoints if self.npoints else 0.0

    # -- processing ---------------------------------------------------------

    def has_fid(self):
        return bool(self.fid)

    def reprocess(self):
        """Rebuild the spectrum: FT (if raw) -> phase -> optional baseline.

        Always starts from the FID or the as-loaded data, so every control in
        the GUI stays non-destructive and re-runnable.
        """
        if self.fid:
            self._transform_fid()
        else:
            self._phase_processed()
        if self.baseline_on:
            self.real = dsp.baseline_correct(self.real, self.baseline_order,
                                             method=self.baseline_method)

    def _transform_fid(self):
        dwell = 1.0 / self.sw_hz
        grpdly = float(self.meta.get("GRPDLY", 0.0) or 0.0)
        fcor = float(self.meta.get("FCOR", 0.5) or 0.5)
        spec = dsp.transform(self.fid, dwell, self.si, lb=self.lb,
                             grpdly=grpdly, fcor=fcor)
        spec = dsp.phase(spec, self.p0, self.p1)
        self.real = [v.real for v in spec]
        self.imag = [v.imag for v in spec]

    def _phase_processed(self):
        """Phase a vendor-processed spectrum, synthesising ``1i`` if absent."""
        if self._base_real is None:
            return
        if self.p0 == 0.0 and self.p1 == 0.0:
            self.real = list(self._base_real)
            self.imag = list(self._base_imag) if self._base_imag else None
            return
        if self._base_imag is None:
            self._base_imag = dsp.hilbert(self._base_real)
        spec = [complex(r, i) for r, i in zip(self._base_real, self._base_imag)]
        spec = dsp.phase(spec, self.p0, self.p1)
        self.real = [v.real for v in spec]
        self.imag = [v.imag for v in spec]

    def complex_spectrum(self):
        """Current spectrum as complex values, for autophase."""
        if self.imag is None:
            self.imag = dsp.hilbert(self.real)
        return [complex(r, i) for r, i in zip(self.real, self.imag)]

    def snapshot_base(self):
        """Remember the as-loaded data so phasing stays non-destructive."""
        self._base_real = list(self.real)
        self._base_imag = list(self.imag) if self.imag else None

    _base_real = None
    _base_imag = None


# ---------------------------------------------------------------------------
# File access that works for both folders and zip archives
# ---------------------------------------------------------------------------


class _Store:
    """Uniform read-only view over a directory or a zip archive."""

    def __init__(self, path):
        self.path = path
        self.zip = None
        if zipfile.is_zipfile(path):
            self.zip = zipfile.ZipFile(path)
            self.names = [n for n in self.zip.namelist() if not n.endswith("/")]
        else:
            self.names = []
            for root, _dirs, files in os.walk(path):
                for f in files:
                    full = os.path.join(root, f)
                    self.names.append(os.path.relpath(full, path).replace(os.sep, "/"))

    def read(self, name):
        if self.zip is not None:
            return self.zip.read(name)
        with open(os.path.join(self.path, name.replace("/", os.sep)), "rb") as fh:
            return fh.read()

    def exists(self, name):
        return name in self.names

    def close(self):
        if self.zip is not None:
            self.zip.close()


# ---------------------------------------------------------------------------
# Bruker
# ---------------------------------------------------------------------------

_PAR_RE = re.compile(r"^##\$([A-Za-z_0-9]+)=\s*(.*)$")


def parse_jcamp_params(text):
    """Parse the ``##$KEY= value`` lines of acqus / procs."""
    params = {}
    for line in text.splitlines():
        m = _PAR_RE.match(line)
        if not m:
            continue
        key, raw = m.group(1), m.group(2).strip()
        if raw.startswith("<") and raw.endswith(">"):
            params[key] = raw[1:-1]
        elif raw.startswith("("):
            params[key] = raw          # array header, kept verbatim
        else:
            try:
                params[key] = float(raw) if ("." in raw or "e" in raw.lower()) else int(raw)
            except ValueError:
                params[key] = raw
    return params


def _unpack(data, dtype, big_endian):
    order = ">" if big_endian else "<"
    if dtype == "int32":
        n = len(data) // 4
        return list(struct.unpack(order + "%di" % n, data[: n * 4]))
    n = len(data) // 8
    return list(struct.unpack(order + "%dd" % n, data[: n * 8]))


def find_bruker_experiments(store):
    """Locate every folder in ``store`` that holds an ``acqus`` file."""
    roots = []
    for name in store.names:
        if name.endswith("acqus") and "/pdata/" not in name:
            roots.append(name[: -len("acqus")])
    return sorted(set(roots))


def available_procnos(store, root):
    """Processing numbers present under an experiment, lowest first.

    TopSpin writes reprocessed data to ``pdata/2``, ``pdata/3`` and so on;
    assuming ``pdata/1`` silently hides that work.
    """
    found = set()
    prefix = root + "pdata/"
    for name in store.names:
        if not name.startswith(prefix):
            continue
        rest = name[len(prefix):].split("/", 1)
        if len(rest) == 2 and rest[0].isdigit():
            found.add(int(rest[0]))
    return sorted(found)


def read_bruker(path, load_fid=True, procno=None):
    """Read all experiments found at ``path`` (folder or zip)."""
    store = _Store(path)
    try:
        specs = []
        for root in find_bruker_experiments(store):
            spec = _read_bruker_one(store, root, path, load_fid, procno)
            if spec is not None:
                specs.append(spec)
        return specs
    finally:
        store.close()


def _read_bruker_one(store, root, path, load_fid, procno=None):
    acqus = parse_jcamp_params(store.read(root + "acqus").decode("latin-1"))

    sfo1 = float(acqus.get("SFO1", 0.0))
    sw_hz = float(acqus.get("SW_h", 0.0))
    if not sfo1 or not sw_hz:
        return None

    # Prefer the processed data; the acquisition parameters alone do not
    # define the referenced ppm axis.  Fall back to whichever processing
    # number actually exists rather than assuming 1.
    candidates = [procno] if procno else available_procnos(store, root) or [1]
    chosen = candidates[0]
    for number in candidates:
        if store.exists("%spdata/%s/1r" % (root, number)):
            chosen = number
            break
    pdata = "%spdata/%s/" % (root, chosen)
    real = imag = None
    sf = sfo1
    offset_ppm = None
    procs = {}
    if store.exists(pdata + "procs") and store.exists(pdata + "1r"):
        procs = parse_jcamp_params(store.read(pdata + "procs").decode("latin-1"))
        big = int(procs.get("BYTORDP", 0)) == 1
        dtype = "float64" if int(procs.get("DTYPP", 0)) == 2 else "int32"
        scale = 2.0 ** float(procs.get("NC_proc", 0))
        real = [v * scale for v in _unpack(store.read(pdata + "1r"), dtype, big)]
        if store.exists(pdata + "1i"):
            imag = [v * scale for v in _unpack(store.read(pdata + "1i"), dtype, big)]
        sf = float(procs.get("SF", sfo1))
        sw_hz = float(procs.get("SW_p", sw_hz))
        offset_ppm = float(procs.get("OFFSET", 0.0))

    fid = None
    if load_fid and store.exists(root + "fid"):
        big = int(acqus.get("BYTORDA", 0)) == 1
        dtype = "float64" if int(acqus.get("DTYPA", 0)) == 2 else "int32"
        raw = _unpack(store.read(root + "fid"), dtype, big)
        td = int(acqus.get("TD", len(raw)))
        raw = raw[: td - (td % 2)]
        fid = [complex(raw[i], raw[i + 1]) for i in range(0, len(raw) - 1, 2)]

    if real is None:
        if fid is None:
            return None
        # No processed data on disk: transform the FID with default settings.
        si = dsp.next_pow2(len(fid))
        spec = dsp.transform(fid, 1.0 / sw_hz, si,
                             lb=float(acqus.get("LB", 0.3) or 0.3),
                             grpdly=float(acqus.get("GRPDLY", 0.0) or 0.0))
        real = [v.real for v in spec]
        imag = [v.imag for v in spec]
        offset_ppm = float(acqus.get("O1", 0.0)) / sf + (sw_hz / sf) / 2.0

    title = ""
    if store.exists(pdata + "title"):
        title = store.read(pdata + "title").decode("latin-1").strip().splitlines()
        title = title[0] if title else ""
    expno = root.rstrip("/").split("/")[-1]
    name = title or "%s [%s]" % (os.path.basename(path), expno)

    meta = {
        "Nucleus": acqus.get("NUC1", ""),
        "Solvent": acqus.get("SOLVENT", ""),
        "Pulse program": acqus.get("PULPROG", ""),
        "Scans": acqus.get("NS", ""),
        "Frequency (MHz)": round(sf, 4),
        "Spectral width (Hz)": round(sw_hz, 2),
        "Temperature (K)": acqus.get("TE", ""),
        "GRPDLY": acqus.get("GRPDLY", 0.0),
        "LB": procs.get("LB", 0.0),
        "FCOR": procs.get("FCOR", 0.5),
        "Format": "Bruker",
        "Experiment": expno,
        "Processing no.": chosen,
    }

    spec = Spectrum(name, real, sf, sw_hz, offset_ppm, imag=imag, fid=fid,
                    meta=meta, source=path)
    spec.p0 = float(procs.get("PHC0", 0.0) or 0.0) * 0  # already applied in 1r
    spec.snapshot_base()
    return spec


# ---------------------------------------------------------------------------
# ACD/Labs .esp
# ---------------------------------------------------------------------------

# Tag numbers observed in ".ESP.( V 1.0 )" files written by ACD/Spectrus.
_ESP_SW_HZ = 0x03
_ESP_OFFSET_HZ = 0x04      # centre of the spectrum, Hz from the reference
_ESP_SF = 0x05             # observe frequency, MHz
_ESP_AQ = 0x06
_ESP_DWELL = 0x07
_ESP_TITLE = 0x0A
_ESP_DATE = 0x0B
_ESP_TEMP = 0x0E
_ESP_SCANS = 0x0F
_ESP_NPOINTS = 0x10        # half the number of stored real points
_ESP_SOLVENT = 0x11
_ESP_PULPROG = 0x19
_ESP_ORIGIN = 0x1B
_ESP_INSTRUMENT = 0x1D
_ESP_OPERATOR = 0x1E


def _esp_header(data):
    """Walk the tag/length/value header and return the tags plus data offset."""
    start = data.find(b"\x03\x04", 0, 512)
    if start < 0:
        raise ValueError("not an ACD .esp spectrum")
    tags = {}
    i = start
    while i < len(data):
        tag = data[i]
        if tag == 0:
            i += 2
            break
        length = data[i + 1]
        tags[tag] = data[i + 2:i + 2 + length]
        i += 2 + length
    return tags, i


def _f32(tags, tag, default=0.0):
    raw = tags.get(tag)
    return struct.unpack("<f", raw)[0] if raw and len(raw) == 4 else default


def _i32(tags, tag, default=0):
    raw = tags.get(tag)
    return struct.unpack("<i", raw)[0] if raw and len(raw) == 4 else default


def _text(tags, tag):
    raw = tags.get(tag)
    return raw.decode("latin-1").strip() if raw else ""


def read_esp(path):
    """Read an ACD/Labs ``.esp`` spectrum.

    The file is an ACD compound document; the spectrum sits immediately after
    the header as two float32 blocks (real then imaginary), each holding twice
    the point count stored in the header, ordered from **low to high ppm**.
    """
    with open(path, "rb") as fh:
        data = fh.read()

    tags, offset = _esp_header(data)
    sw_hz = _f32(tags, _ESP_SW_HZ)
    sf = _f32(tags, _ESP_SF)
    centre_hz = _f32(tags, _ESP_OFFSET_HZ)
    npts = _i32(tags, _ESP_NPOINTS) * 2
    if not (sf and sw_hz and npts):
        raise ValueError("incomplete .esp header")

    need = offset + npts * 4
    if need > len(data):
        raise ValueError("truncated .esp spectrum")
    real = list(struct.unpack_from("<%df" % npts, data, offset))
    imag = None
    if offset + 2 * npts * 4 <= len(data):
        imag = list(struct.unpack_from("<%df" % npts, data, offset + npts * 4))

    # Stored low-to-high ppm; flip to the high-ppm-first convention.
    real.reverse()
    if imag:
        imag.reverse()

    sw_ppm = sw_hz / sf
    right_ppm = centre_hz / sf - sw_ppm / 2.0
    offset_ppm = right_ppm + (npts - 1) * sw_ppm / npts

    meta = {
        "Solvent": _text(tags, _ESP_SOLVENT),
        "Pulse program": _text(tags, _ESP_PULPROG),
        "Scans": struct.unpack("<h", tags[_ESP_SCANS])[0] if _ESP_SCANS in tags else "",
        "Frequency (MHz)": round(sf, 4),
        "Spectral width (Hz)": round(sw_hz, 2),
        "Temperature (C)": round(_f32(tags, _ESP_TEMP), 2),
        "Date": _text(tags, _ESP_DATE),
        "Instrument": _text(tags, _ESP_INSTRUMENT),
        "Operator": _text(tags, _ESP_OPERATOR),
        "Origin": _text(tags, _ESP_ORIGIN),
        "Format": "ACD .esp",
    }
    name = _text(tags, _ESP_TITLE) or os.path.splitext(os.path.basename(path))[0]

    spec = Spectrum(name.strip(), real, sf, sw_hz, offset_ppm, imag=imag,
                    meta=meta, source=path)
    spec.snapshot_base()
    return [spec]


# ---------------------------------------------------------------------------
# JCAMP-DX
# ---------------------------------------------------------------------------

# ASDF compression alphabets (JCAMP-DX 4.24 / 5.01).
_SQZ = {"@": 0, "A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6, "G": 7,
        "H": 8, "I": 9, "a": -1, "b": -2, "c": -3, "d": -4, "e": -5,
        "f": -6, "g": -7, "h": -8, "i": -9}
_DIF = {"%": 0, "J": 1, "K": 2, "L": 3, "M": 4, "N": 5, "O": 6, "P": 7,
        "Q": 8, "R": 9, "j": -1, "k": -2, "l": -3, "m": -4, "n": -5,
        "o": -6, "p": -7, "q": -8, "r": -9}
_DUP = {"S": 1, "T": 2, "U": 3, "V": 4, "W": 5, "X": 6, "Y": 7, "Z": 8, "s": 9}


def _decode_asdf(line):
    """Decode one ASDF ordinate line into ``(values, ended_in_dif)``.

    Handles all of PAC, SQZ, DIF and DUP, including decimal points, which
    matter because in SQZ and DIF files the abscissa runs straight into the
    first ordinate with no separator (``5000.03B1399T...``).
    """
    values = []
    i = 0
    n = len(line)
    last_was_dif = False
    while i < n:
        ch = line[i]
        if ch in " ,\t":
            i += 1
            continue

        if ch in _DUP:
            count = _DUP[ch]
            i += 1
            while i < n and line[i].isdigit():
                count = count * 10 + int(line[i])
                i += 1
            if not values:
                continue
            # In DIF form a DUP repeats the last difference, otherwise the
            # last value itself.
            step = (values[-1] - values[-2]
                    if last_was_dif and len(values) > 1 else 0.0)
            for _ in range(count - 1):
                values.append(values[-1] + step)
            continue

        if ch in _DIF or ch in _SQZ or ch in "+-" or ch.isdigit() or ch == ".":
            is_dif = False
            lead = ""
            negative = False
            if ch in _DIF:
                lead, is_dif, negative = str(abs(_DIF[ch])), True, _DIF[ch] < 0
                i += 1
            elif ch in _SQZ:
                lead, negative = str(abs(_SQZ[ch])), _SQZ[ch] < 0
                i += 1
            else:
                if ch in "+-":
                    negative = ch == "-"
                    i += 1
            digits = ""
            while i < n and (line[i].isdigit() or line[i] == "."):
                digits += line[i]
                i += 1
            # No exponent handling on purpose: "E" is the SQZ digit 5, so
            # reading "1946E434" as 1946e434 would silently produce infinity.
            # JCAMP data tables never use scientific notation.
            try:
                num = float(lead + digits) if (lead or digits) else 0.0
            except ValueError:
                num = 0.0
            if negative:
                num = -num
            if is_dif and values:
                values.append(values[-1] + num)
            else:
                values.append(num)
            last_was_dif = is_dif
            continue
        i += 1
    return values, last_was_dif


def _jcamp_key(raw):
    """Normalise a label so ``##.OBSERVE FREQUENCY`` and ``##OBSERVEFREQUENCY``
    end up as the same key."""
    return "".join(ch for ch in raw.upper() if ch.isalnum())


def _decode_table(lines, start):
    """Read the ordinate rows that follow a data-table label.

    The first number on each row is the abscissa and is discarded; the axis is
    reconstructed from FIRSTX/LASTX instead, which is both simpler and what
    the format guarantees.
    """
    ydata = []
    previous_was_dif = False
    for raw in lines[start:]:
        if raw.startswith("##"):
            break
        body = raw.split("$$")[0].strip()
        if not body:
            continue
        values, ended_in_dif = _decode_asdf(body)
        if len(values) < 2:
            continue
        chunk = values[1:]
        # A DIF row restates the previous row's last ordinate as a checkpoint.
        if previous_was_dif and ydata and abs(chunk[0] - ydata[-1]) < 1e-6:
            chunk = chunk[1:]
        ydata.extend(chunk)
        previous_was_dif = ended_in_dif
    return ydata


class _JcampBlock:
    def __init__(self):
        self.header = {}
        self.tables = []        # (label_value, first_data_line_index)

    def get(self, key, default=None):
        return self.header.get(_jcamp_key(key), default)

    def number(self, key, default=0.0):
        try:
            return float(str(self.get(key, default)).split()[0])
        except (TypeError, ValueError, IndexError):
            return default


def _split_blocks(lines):
    """Split a JCAMP file into blocks; LINK files hold several."""
    blocks = []
    current = None
    for index, raw in enumerate(lines):
        if not raw.startswith("##"):
            continue
        label, _, value = raw[2:].partition("=")
        key = _jcamp_key(label)
        value = value.split("$$")[0].strip()
        if key == "TITLE":
            current = _JcampBlock()
            blocks.append(current)
        if current is None:
            current = _JcampBlock()
            blocks.append(current)
        current.header[key] = value
        if key in ("XYDATA", "DATATABLE", "XYPOINTS"):
            current.tables.append((value, index + 1))
    return blocks


def _pick_spectrum_block(blocks):
    """Choose the block that actually holds spectral data."""
    scored = []
    for block in blocks:
        if not block.tables:
            continue
        kind = str(block.get("DATATYPE", "")).upper()
        npoints = block.number("NPOINTS", 0)
        score = npoints
        if "SPECTRUM" in kind:
            score += 1e9
        if "FID" in kind:
            score += 5e8
        if "ASSIGNMENT" in kind or "PEAK" in kind:
            score -= 1e9
        scored.append((score, block))
    if not scored:
        return None
    return max(scored, key=lambda item: item[0])[1]


def _ntuples_factors(block):
    """Map NTUPLES symbols to their per-symbol FACTOR / FIRST / LAST."""
    def split_list(key):
        raw = block.get(key, "")
        return [part.strip() for part in str(raw).split(",")]

    symbols = split_list("SYMBOL")
    out = {}
    for name in ("FACTOR", "FIRST", "LAST", "VARDIM"):
        values = split_list(name)
        for position, symbol in enumerate(symbols):
            if position < len(values):
                try:
                    out.setdefault(symbol, {})[name] = float(values[position])
                except ValueError:
                    pass
    return out


def read_jcamp(path):
    """Read a JCAMP-DX spectrum or FID.

    Understands single-block files, multi-block LINK files (the spectrum block
    is picked automatically) and NTUPLES data tables with separate real and
    imaginary pages.
    """
    with open(path, "r", encoding="latin-1") as fh:
        lines = fh.read().splitlines()

    blocks = _split_blocks(lines)
    block = _pick_spectrum_block(blocks)
    if block is None:
        raise ValueError("no spectral data table in this JCAMP file")

    sf = block.number("OBSERVEFREQUENCY", 0.0) or block.number("BF1", 0.0)
    yfactor = block.number("YFACTOR", 1.0) or 1.0
    xunits = str(block.get("XUNITS", "")).upper()
    kind = str(block.get("DATATYPE", "")).upper()

    # -- gather the data pages
    pages = {}
    for label, start in block.tables:
        symbol = "Y"
        if "(R.." in label.replace(" ", ""):
            symbol = "R"
        elif "(I.." in label.replace(" ", ""):
            symbol = "I"
        pages[symbol] = _decode_table(lines, start)

    ntuples = _ntuples_factors(block) if block.get("NTUPLES") else {}
    if ntuples:
        for symbol, values in pages.items():
            factor = ntuples.get(symbol, {}).get("FACTOR", 1.0)
            pages[symbol] = [v * factor for v in values]
        firstx = ntuples.get("X", {}).get("FIRST", block.number("FIRSTX", 0.0))
        lastx = ntuples.get("X", {}).get("LAST", block.number("LASTX", 0.0))
    else:
        for symbol, values in pages.items():
            pages[symbol] = [v * yfactor for v in values]
        firstx = block.number("FIRSTX", 0.0)
        lastx = block.number("LASTX", 0.0)

    real = pages.get("R") or pages.get("Y") or []
    imag = pages.get("I")
    if not real:
        raise ValueError("the JCAMP data table decoded to nothing")

    meta = {
        "Solvent": block.get("SOLVENTNAME", ""),
        "Nucleus": block.get("OBSERVENUCLEUS", ""),
        "Pulse program": block.get("PULSESEQUENCE", ""),
        "Origin": block.get("ORIGIN", ""),
        "Format": "JCAMP-DX %s" % block.get("JCAMPDX", ""),
        "X units": xunits or "?",
    }
    name = block.get("TITLE") or os.path.basename(path)

    # -- a FID is transformed rather than displayed as-is
    if "FID" in kind and imag and sf:
        fid = [complex(r, i) for r, i in zip(real, imag)]
        duration = abs(lastx - firstx)
        sw_hz = (len(fid) - 1) / duration if duration else 0.0
        if sw_hz:
            si = dsp.next_pow2(len(fid))
            spec = dsp.transform(fid, 1.0 / sw_hz, si, lb=0.3)
            offset_ppm = (sw_hz / sf) / 2.0
            meta["Format"] += " (FID, transformed)"
            result = Spectrum(name, [v.real for v in spec], sf, sw_hz,
                              offset_ppm, imag=[v.imag for v in spec],
                              fid=fid, meta=meta, source=path)
            result.snapshot_base()
            return [result]

    # -- frequency-domain data
    if firstx < lastx:
        real.reverse()
        if imag:
            imag.reverse()
        firstx, lastx = lastx, firstx

    if xunits.startswith("HZ") and sf:
        firstx, lastx = firstx / sf, lastx / sf     # Hz -> ppm
    elif not xunits.startswith("PPM"):
        # Not an NMR abscissa (1/CM, NANOMETERS, ...); show the native units.
        sf = sf or 1.0

    n = len(real)
    span = abs(firstx - lastx)
    # FIRSTX/LASTX mark the first and last points, so the step is span/(n-1);
    # Spectrum stores a width whose step is span/n.
    sw_ppm = span * n / (n - 1) if n > 1 else span
    if not sf:
        sf = 1.0

    spec = Spectrum(name, real, sf, sw_ppm * sf, firstx, imag=imag,
                    meta=meta, source=path)
    spec.snapshot_base()
    return [spec]


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def load(path):
    """Load any supported path and return a list of :class:`Spectrum`."""
    lower = path.lower()
    if lower.endswith(".esp"):
        specs = read_esp(path)
    elif lower.endswith((".jdx", ".dx", ".jcamp")):
        specs = read_jcamp(path)
    elif os.path.isdir(path) or zipfile.is_zipfile(path):
        specs = read_bruker(path)
        if not specs:
            raise ValueError("no Bruker experiment (acqus) found in %s"
                             % os.path.basename(path))
    else:
        raise ValueError("unsupported file type: %s" % os.path.basename(path))

    # Sessions refer back to a spectrum by its position within its source.
    for index, spec in enumerate(specs):
        spec.source_index = index
    return specs
