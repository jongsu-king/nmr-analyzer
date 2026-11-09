#!/usr/bin/env bash
# Package the one-folder build into a .tar.gz (with install.sh) and a .deb.
# Usage: packaging/linux/build_linux.sh 1.6.0   (after pyinstaller, from the repository root)
set -euo pipefail
VERSION="$1"
mkdir -p out
ID=nmr-analyzer
SLUG=NMR-Analyzer
BUILD=NMR-Analyzer           # PyInstaller output folder and executable name
NAME="NMR Analyzer"
SUMMARY="Viewer and analyser for 1D and 2D NMR spectra"

# --- portable archive with a per-user installer
STAGE=$(mktemp -d)/$SLUG
cp -R "dist-app/$BUILD" "$STAGE"
for f in install.sh uninstall.sh; do
  sed -e "s|@ID@|$ID|g" -e "s|@EXE@|$BUILD|g" "packaging/linux/$f" > "$STAGE/$f"
done
cp packaging/icon.png "$STAGE/$ID.png"
sed -e "s|@NAME@|$NAME|" -e "s|@SUMMARY@|$SUMMARY|" -e "s|@ID@|$ID|" packaging/linux/app.desktop > "$STAGE/$ID.desktop"
cp LICENSE README.md "$STAGE/"
chmod +x "$STAGE/install.sh" "$STAGE/uninstall.sh"
tar -C "$(dirname "$STAGE")" -czf "out/${SLUG}-${VERSION}-linux-x86_64.tar.gz" "$SLUG"

# --- Debian/Ubuntu package
PKG=$(mktemp -d)/${ID}_${VERSION}_amd64
mkdir -p "$PKG/DEBIAN" "$PKG/opt" "$PKG/usr/bin" "$PKG/usr/share/applications" \
         "$PKG/usr/share/icons/hicolor/256x256/apps" "$PKG/usr/share/doc/$ID"
cp -R "dist-app/$BUILD" "$PKG/opt/$ID"
ln -s "/opt/$ID/$BUILD" "$PKG/usr/bin/$ID"
sed -e "s|@NAME@|$NAME|" -e "s|@SUMMARY@|$SUMMARY|" -e "s|@ID@|$ID|" -e "s|@EXEC@|/opt/$ID/$BUILD|" \
    packaging/linux/app.desktop > "$PKG/usr/share/applications/$ID.desktop"
cp packaging/icon.png "$PKG/usr/share/icons/hicolor/256x256/apps/$ID.png"
cp LICENSE "$PKG/usr/share/doc/$ID/copyright"
cat > "$PKG/DEBIAN/control" <<CONTROL
Package: $ID
Version: $VERSION
Section: science
Priority: optional
Architecture: amd64
Installed-Size: $(du -sk "$PKG" | cut -f1)
Depends: libc6 (>= 2.35), libx11-6, libxft2, libxss1, libfontconfig1
Maintainer: Jongsu Lim <313024960+jongsu-king@users.noreply.github.com>
Homepage: https://github.com/jongsu-king/nmr-analyzer
Description: $SUMMARY
 Reads Bruker, ACD/Labs .esp and JCAMP-DX data; processing, peak picking,
 integration, line-shape fitting, 2D contour display and a structure check.
CONTROL
dpkg-deb --build --root-owner-group "$PKG" "out/${ID}_${VERSION}_amd64.deb"
ls -la out
