; Inno Setup script for the Windows installer.
;
; The payload is the single PyInstaller executable, so there is nothing to
; unpack beyond copying it into Program Files and creating the shortcuts.
; AppVersion is passed in by the build workflow with /DAppVersion=...

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

#define AppName "NMR Analyzer"
#define AppExe "NMR-Analyzer.exe"
#define Publisher "jongsu-king"
#define AppUrl "https://github.com/jongsu-king/nmr-analyzer"

[Setup]
AppId={{7B3C1F42-9A5E-4C2D-8E71-2F0A6D4B9C13}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#Publisher}
AppPublisherURL={#AppUrl}
AppSupportURL={#AppUrl}/issues
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
OutputBaseFilename=NMR-Analyzer-Setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
SetupIconFile=icon.ico
UninstallDisplayIcon={app}\{#AppExe}
; Per-user install needs no administrator rights.
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesInstallIn64BitMode=x64compatible

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "korean"; MessagesFile: "compiler:Languages\Korean.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; \
  GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "..\dist-app\{#AppExe}"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\README.md"; DestDir: "{app}"; DestName: "README.md"; \
  Flags: ignoreversion
Source: "..\LICENSE"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{group}\{cm:UninstallProgram,{#AppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; \
  Tasks: desktopicon

[Registry]
; Offer the app in "Open with" for the formats it reads.
Root: HKA; Subkey: "Software\Classes\.esp\OpenWithProgids"; \
  ValueType: string; ValueName: "NMRAnalyzer.spectrum"; ValueData: ""; \
  Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\.jdx\OpenWithProgids"; \
  ValueType: string; ValueName: "NMRAnalyzer.spectrum"; ValueData: ""; \
  Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\.nmrs\OpenWithProgids"; \
  ValueType: string; ValueName: "NMRAnalyzer.session"; ValueData: ""; \
  Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\NMRAnalyzer.spectrum"; \
  ValueType: string; ValueName: ""; ValueData: "NMR spectrum"; \
  Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\NMRAnalyzer.session"; \
  ValueType: string; ValueName: ""; ValueData: "NMR Analyzer session"; \
  Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\NMRAnalyzer.spectrum\shell\open\command"; \
  ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"" ""%1"""
Root: HKA; Subkey: "Software\Classes\NMRAnalyzer.session\shell\open\command"; \
  ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"" ""%1"""

[Run]
Filename: "{app}\{#AppExe}"; \
  Description: "{cm:LaunchProgram,{#StringChange(AppName, '&', '&&')}}"; \
  Flags: nowait postinstall skipifsilent
