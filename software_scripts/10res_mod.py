

with open("data/9mirnas/finalTable.txt") as res_tab:
    for eli in res_tab.readlines():
        eli = eli.strip().split("\t")
        if eli[21].startswith("zma"):
