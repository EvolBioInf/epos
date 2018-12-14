git describe | sed -r 's/^[vV]//; s/-g(.......)/\\ \\(\1\\)\\/' | tr -d '\n'
