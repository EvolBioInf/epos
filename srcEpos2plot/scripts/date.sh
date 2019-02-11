git log       |
    grep Date |
    head -n 1 |
    sed -E 's/Date: +//g;s/ \+.+//;s/ /_/g'
