# -*- mode: python -*-

block_cipher = None


a = Analysis(['changes_to_pine_final.py'],
             pathex=['C:\\Users\\GoJ1\\Documents\\cytoscape1'],
             binaries=[],
             datas=[],
             hiddenimports=['igraph.vendor.texttable'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='changes_to_pine_final',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='changes_to_pine_final')
