# -*- mode: python -*-

block_cipher = None


a = Analysis(['pyfdap_app.py'],
             pathex=['/home/alex_loc/Documents/Research/PyFDAP/pyfdap'],
             binaries=None,
             datas=None,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='pyfdap_Linux_3.13.0-61-generic_64bit',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='logo/pyfdap_icon.ico')
