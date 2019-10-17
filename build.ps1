set "PATH=%PATH%;C:\Windows\System32\downlevel;"
env\Scripts\activate.ps1
pyinstaller --clean --hidden-import igraph.vendor.texttable .\pine_2.py