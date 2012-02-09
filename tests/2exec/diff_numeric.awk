# Processes the diff of two files and ignores the small errors in numbers -
# either small absolute or relative errors - controlled by 'abs_tol' and 'rel_tol'.
# These values can be either specified outside by e.g. '-v abs_tol=1e-10' or assume
# default values specified below.
# Numbers should be separated from other text by spaces, '=', ',', or brackets -
# controlled by 'FS'. Output happens only if differences are found.

BEGIN {
	i1=0
	i2=0
	number="^[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?$"
	FS="[[:space:]=,{}()[\\]]"
	if (abs_tol !~ number) abs_tol=1e-15
	if (rel_tol !~ number) rel_tol=1e-8
}

function abs(a) {
	return a>0 ? a : -a
}

function relerr(a,b,  c) {
	c=(abs(a)+abs(b))/2;
	return (c==0) ? 0 : abs(a-b)/c
}

/^[0-9]/ {
	if (i1>i2) {
		printf "First file contains unmatched line:\n%s\n", ref[i2+1]
		exit 3
	}
	i1=0
	i2=0
}

/^< / {
	i1++
	ref[i1]=substr($0,3)
}

/^> / {
	i2++
	test=substr($0,3)
	if (i2>i1) {
		printf "Second file contains unmatched line:\n%s\n", test
		exit 3
	}
	if (split(ref[i2],p1) != split(test,p2)) {
		printf "Different number of elements in lines:\n%s\nand\n%s\n", ref[i2], test
		exit 2
	}
	for (j in p1) {
		if (p1[j]!=p2[j]) {
			if (p1[j] !~ number || p2[j] !~ number) {
				printf "Difference between elements: '%s' and '%s'\n", p1[j], p2[j]
				exit 2
			}
			else if (abs(p1[j]-p2[j])>abs_tol && relerr(p1[j],p2[j])>rel_tol) {
				printf "Significant difference between numbers: '%s' and '%s' (abs_tol='%s', rel_tol='%s')\n", p1[j], p2[j], abs_tol, rel_tol
				exit 1
			}
		}
	}
}