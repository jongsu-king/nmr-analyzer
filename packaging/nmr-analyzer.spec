# PyInstaller spec, shared by the Windows, macOS and Linux builds.
#
# The application is pure standard library, so there is nothing to collect
# beyond tkinter itself.  Everything scientific is excluded explicitly: a
# stray numpy or matplotlib on the build machine would otherwise be pulled in
# and quadruple the download for no benefit.

import os
import re
import sys

NAME = "NMR-Analyzer"
ICON_DIR = SPECPATH
# PyInstaller resolves paths in a spec relative to the spec file, so derive
# the repository root from SPECPATH rather than from the working directory.
ROOT = os.path.abspath(os.path.join(SPECPATH, os.pardir))

EXCLUDES = [
    "numpy", "scipy", "matplotlib", "pandas", "PIL", "IPython",
    "pytest", "setuptools", "pip", "wheel", "distutils",
    "sqlite3", "xmlrpc", "pydoc_data", "unittest",
]
# (http, email and ssl stay: the automatic updater talks HTTPS to GitHub.)
with open(os.path.join(ROOT, "nmranalyzer", "__init__.py"), encoding="utf-8") as fh:
    VERSION = re.search(r'__version__ = "([^"]+)"', fh.read()).group(1)

analysis = Analysis(
    ["entry.py"],
    pathex=[ROOT],
    binaries=[],
    datas=[],
    hiddenimports=["tkinter", "tkinter.ttk", "tkinter.filedialog",
                   "tkinter.messagebox", "tkinter.colorchooser",
                   "nmranalyzer.selfupdate"],
    hookspath=[],
    runtime_hooks=[],
    excludes=EXCLUDES,
    noarchive=False,
)
pyz = PYZ(analysis.pure)

# One folder on every platform: it starts faster than a self-unpacking single
# file, and the in-app updater replaces the folder (or the .app) as a whole.
exe = EXE(
    pyz, analysis.scripts, [],
    exclude_binaries=True,
    name=NAME,
    console=False,
    icon=os.path.join(ICON_DIR, "icon.icns" if sys.platform == "darwin" else "icon.ico"),
)
collect = COLLECT(
    exe, analysis.binaries, analysis.datas,
    strip=False, upx=False, name=NAME,
)
if sys.platform == "darwin":
    app = BUNDLE(
        collect,
        name="%s.app" % NAME,
        icon=os.path.join(ICON_DIR, "icon.icns"),
        bundle_identifier="io.github.jongsuking.nmranalyzer",
        info_plist={
            "CFBundleName": "NMR Analyzer",
            "CFBundleDisplayName": "NMR Analyzer",
            "CFBundleShortVersionString": VERSION,
            "CFBundleVersion": VERSION,
            "NSHighResolutionCapable": True,
            "LSMinimumSystemVersion": "11.0",
            "CFBundleDocumentTypes": [{
                "CFBundleTypeName": "NMR data",
                "CFBundleTypeRole": "Viewer",
                "LSItemContentTypes": ["public.data"],
                "CFBundleTypeExtensions": ["esp", "jdx", "dx", "nmrs", "zip"],
            }],
        },
    )
