function [ppos61, flag_data_min] = getOptimumPpos61()

ppos61 = [1 11 91 104 191 236 262 269 275 381 393 478 519 532 583 597 607 654 733 829 872 899 934 980 1041 1061 1202 1211 1245 1259 1357 1371 1470 1485 1510 1516 1521 1526 1551 1576 1588 1679 1686 1710 1732 1743 1746 1768 1777 1845 1852 1855 1856 1872 1884 1929 1938 1942 1987 1995 2014;44 50 76 100 204 294 317 318 353 363 415 437 446 449 484 501 641 661 663 701 708 715 726 772 785 828 858 882 884 910 920 986 1050 1110 1166 1168 1228 1291 1352 1367 1374 1453 1486 1592 1597 1674 1678 1725 1764 1766 1781 1812 1849 1858 1879 1892 1896 1947 2034 2039 2043;7 136 166 175 198 247 250 310 406 436 441 447 477 482 483 575 585 615 657 691 705 740 744 762 774 838 849 870 954 978 987 996 1004 1027 1038 1106 1183 1187 1216 1231 1275 1287 1305 1328 1334 1393 1431 1505 1522 1587 1657 1718 1719 1769 1828 1835 1900 1996 1999 2016 2020;37 165 218 237 347 350 391 397 463 507 512 599 617 643 644 648 669 728 734 781 795 836 861 924 932 941 959 960 1037 1074 1125 1161 1169 1175 1188 1214 1230 1263 1269 1283 1308 1419 1451 1504 1534 1538 1590 1745 1756 1759 1792 1803 1843 1859 1883 1921 1923 1934 1976 1993 2037;2 35 43 55 59 67 112 134 154 203 362 452 556 642 647 670 710 713 763 779 796 797 822 835 837 862 890 894 903 944 949 1059 1093 1097 1124 1134 1163 1174 1179 1197 1237 1248 1257 1300 1339 1353 1399 1426 1434 1458 1467 1507 1575 1631 1730 1758 1873 1903 1936 2022 2023;3 26 29 56 89 150 208 221 224 226 261 338 373 377 414 558 564 565 582 584 799 898 935 1066 1078 1082 1096 1116 1122 1126 1141 1151 1191 1201 1262 1361 1372 1389 1469 1476 1491 1492 1506 1550 1566 1593 1595 1636 1642 1748 1820 1829 1842 1870 1899 1910 1969 1981 1989 2001 2009;95 124 125 135 177 209 268 355 360 389 399 461 520 567 624 633 732 790 847 897 921 931 963 1011 1060 1075 1111 1150 1200 1203 1224 1251 1252 1293 1309 1324 1382 1477 1483 1495 1548 1552 1567 1612 1669 1671 1692 1760 1838 1874 1906 1908 1935 1944 1950 1967 1973 1977 1985 2013 2029;63 129 186 211 259 263 286 293 351 368 396 431 495 526 529 530 544 577 610 682 926 977 1047 1067 1077 1083 1100 1143 1171 1223 1249 1253 1299 1306 1318 1366 1405 1415 1436 1440 1462 1475 1498 1559 1572 1573 1606 1611 1666 1682 1737 1795 1822 1827 1901 1918 1924 1956 1979 2002 2017;64 101 105 114 131 149 163 187 219 251 267 271 282 357 421 443 469 540 542 570 572 586 596 687 760 764 805 820 852 868 871 925 958 997 1002 1054 1058 1068 1085 1109 1270 1271 1354 1474 1489 1546 1547 1549 1602 1607 1717 1787 1808 1888 1893 1932 1957 1958 1978 2019 2044;9 15 19 51 132 133 144 289 290 292 321 343 416 423 487 489 562 631 665 722 730 731 746 759 776 783 823 876 878 886 889 911 993 1005 1013 1018 1070 1112 1140 1167 1236 1265 1280 1319 1402 1407 1461 1528 1577 1617 1649 1697 1753 1790 1797 1802 1851 1871 1881 1962 2018];

load('flag_dataK10_Np61.mat');
