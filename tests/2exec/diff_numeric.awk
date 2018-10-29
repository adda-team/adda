# Processes the diff of two files and ignores the small errors in numbers - either small absolute or relative errors -
# controlled by 'abs_tol' and 'rel_tol'. These values can be either specified outside by e.g. '-v abs_tol=1e-10' or
# assume default values specified below. Numbers should be separated from other text by spaces, '=', ',', ":" or 
# brackets - controlled by 'FS'. Complex numbers are also supported, for them real and imaginary parts are processed
# independently. Output happens only if differences are found.

BEGIN {
	i1=0
	i2=0
	num="[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?"
	number="^" num "$"
	complex="^" num num "[iI]$"
	# FS="[[:space:]=,{}()[\\]]"
	# A more robust assignment of FS is used to be compatible with mawk
	FS="[ \\t=,{}\\(\\)[\\]:]"
	if (abs_tol !~ number) abs_tol=1e-15
	if (rel_tol !~ number) rel_tol=1e-8
	exit_invoked=0
}

function abs(a) {
	return a>0 ? a : -a
}

function relerr(a,b,  c) {
	c=(abs(a)+abs(b))/2;
	return (c==0) ? 0 : abs(a-b)/c
}

function stop(a) {
	exit_invoked=1
	exit a
}

function parseComplex(str,c) {
	match(str,num);
	c["re"]=substr(str,RSTART+1,RLENGTH-1)
	match(str,num "[iI]$");
	c["im"]=substr(str,RSTART+1,RLENGTH-2)
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
		stop(3)
	}
	if (split(ref[i2],p1) != split(test,p2)) {
		printf "Different number of elements in lines:\n%s\nand\n%s\n", ref[i2], test
		stop(2)
	}
	for (j in p1) {
		if (p1[j]!=p2[j]) {
			if (p1[j] ~ number && p2[j] ~ number) {
				if (abs(p1[j]-p2[j])>abs_tol && relerr(p1[j],p2[j])>rel_tol) {
					# wrapping long lines produces some errors on Windows, so keep it like this
					printf "Significant difference between numbers: '%s' and '%s' (abs_tol='%s', rel_tol='%s')\n", p1[j], p2[j], abs_tol, rel_tol
					stop(1)
				}
			}
			else if (p1[j] ~ complex && p2[j] ~ complex) {
				parseComplex(p1[j],c1)
				parseComplex(p2[j],c2)
				if ( (abs(c1["re"]-c2["re"])>abs_tol && relerr(c1["re"],c2["re"])>rel_tol) || (abs(c1["im"]-c2["im"])>abs_tol && relerr(c1["im"],c2["im"])>rel_tol) ) {
					printf "Significant difference between complex numbers: '%s' and '%s' (abs_tol='%s', rel_tol='%s')\n", p1[j], p2[j], abs_tol, rel_tol
					stop(1)
				}
			}
			else {
				printf "Difference between elements: '%s' and '%s'\n", p1[j], p2[j]
				stop(2)
			}
		}
	}
}

/^[0-9]/ {
	if (i1>i2) {
		printf "First file contains unmatched line:\n%s\n", ref[i2+1]
		stop(3)
	}
	i1=0
	i2=0
}

END {
	if (!exit_invoked) {
		if (i1>i2) {
			printf "First file contains unmatched line:\n%s\n", ref[i2+1]
			exit 3
		}
	}
}