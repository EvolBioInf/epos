git describe | sed -r 's/^[vV]//; s/(.)$/\1\\/; s/-g(.......)/\\ \\(\1\\)/' | tr -d '\n'
