set "PATH=%PATH%;C:\Windows\System32\downlevel;"
pyinstaller --hidden-import igraph.vendor.texttable .\pine_2.py