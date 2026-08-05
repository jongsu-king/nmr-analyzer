# PyInstaller spec, shared by the Windows, macOS and Linux builds.
#
# The application is pure standard library, so there is nothing to collect
# beyond tkinter itself.  Everything scientific is excluded explicitly: a
# stray numpy or matplotlib on the build machine would otherwise be pulled in
# and quadruple the download for no benefit.

import os
import sys

NAME = "NMR-Analyzer"
ICON_DIR = SPECPATH
# PyInstaller resolves paths in a spec relative to the spec file, so derive
# the repository root from SPECPATH rather than from the working directory.
ROOT = os.path.abspath(os.path.join(SPECPATH, os.pardir))

EXCLUDES = [
    "numpy", "scipy", "matplotlib", "pandas", "PIL", "IPython",
    "pytest", "setuptools", "pip", "wheel", "distutils",
    "sqlite3", "email", "http", "xmlrpc", "pydoc_data", "unittest",
]

analysis = Analysis(
    ["entry.py"],
    pathex=[ROOT],
    binaries=[],
    datas=[],
    hiddenimports=["tkinter", "tkinter.ttk", "tkinter.filedialog",
                   "tkinter.messagebox", "tkinter.colorchooser"],
    hookspath=[],
    runtime_hooks=[],
    excludes=EXCLUDES,
    noarchive=False,
)
pyz = PYZ(analysis.pure)

if sys.platform == "darwin":
    # One directory, wrapped in a .app -- macOS refuses to notarise, and
    # Gatekeeper is unhappy with, a single-file binary that unpacks itself.
    exe = EXE(
        pyz, analysis.scripts, [],
        exclude_binaries=True,
        name=NAME,
        console=False,
        icon=os.path.join(ICON_DIR, "icon.icns"),
    )
    collect = COLLECT(
        exe, analysis.binaries, analysis.datas,
        strip=False, upx=False, name=NAME,
    )
    app = BUNDLE(
        collect,
        name="%s.app" % NAME,
        icon=os.path.join(ICON_DIR, "icon.icns"),
        bundle_identifier="io.github.jongsuking.nmranalyzer",
        info_plist={
            "CFBundleName": "NMR Analyzer",
            "CFBundleDisplayName": "NMR Analyzer",
            "CFBundleShortVersionString": os.environ.get("NMR_VERSION", "0.0.0"),
            "CFBundleVersion": os.environ.get("NMR_VERSION", "0.0.0"),
            "NSHighResolutionCapable": True,
            "LSMinimumSystemVersion": "10.13",
            "CFBundleDocumentTypes": [{
                "CFBundleTypeName": "NMR data",
                "CFBundleTypeRole": "Viewer",
                "LSItemContentTypes": ["public.data"],
                "CFBundleTypeExtensions": ["esp", "jdx", "dx", "nmrs", "zip"],
            }],
        },
    )
else:
    # A single self-contained file is the friendliest thing to hand someone
    # on Windows and Linux.
    exe = EXE(
        pyz, analysis.scripts, analysis.binaries, analysis.datas, [],
        name=NAME,
        console=False,
        strip=False,
        upx=False,
        icon=os.path.join(ICON_DIR, "icon.ico"),
    )
