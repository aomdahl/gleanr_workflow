#!/bin/awk -f

# Set the field separator to any whitespace
BEGIN {
    FS = "[ \t]+"
}

{

	maf = ((1-($10 + 0.0))<$10) ? (1-($10 + 0.0)): $10

	if ((NR == 1) || (($7 == "true") && ($8 + 0.0 > 0.9) && (maf > 0.01) && (length($3)==1) && (length($4)==1))) {
		if (NR == 1) {
			print $0"\tMAF"
		}
		if (NR != 1) {
			print $0"\t"maf
		}
	}
}
