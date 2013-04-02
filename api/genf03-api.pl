#!/usr/bin/perl -w
# Generate Fortran 2003 interfaces from a sequence of C function declarations
# of the form (one per line):
#     extern <type> <name>(...args...)
#     extern <type> <name>(...args...)
#     ...
# with no line breaks within a given function.  (It's too much work to
# write a general parser, since we just have to handle FFTW's header files.)

sub canonicalize_type {
    my($type);
    ($type) = @_;
    $type =~ s/ +/ /g;
    $type =~ s/^ //;
    $type =~ s/ $//;
    $type =~ s/([^\* ])\*/$1 \*/g;
    return $type;
}

sub canonicalize_name {
    my($name);
    ($name) = @_;
    # Since Fortran isn't case sensitive we must distiguish between N and n
    if ($name eq "n") {
      $name = "Nos";
    }
    return $name;
}

# C->Fortran map of supported return types
%return_types = (
    "int" => "integer(C_INT)",
    "unsigned" => "integer(C_INT)",
    "ptrdiff_t" => "integer(C_INTPTR_T)",
    "size_t" => "integer(C_SIZE_T)",
    "double" => "real(C_DOUBLE)",
    "float" => "real(C_FLOAT)",
    "long double" => "real(C_LONG_DOUBLE)",
    "float128__" => "real(16)",
    "pfft_plan" => "type(C_PTR)",
    "pfftf_plan" => "type(C_PTR)",
    "pfftl_plan" => "type(C_PTR)",
    "pfftq_plan" => "type(C_PTR)",
    "pfft_gcplan" => "type(C_PTR)",
    "pfftf_gcplan" => "type(C_PTR)",
    "pfftl_gcplan" => "type(C_PTR)",
    "pfftq_gcplan" => "type(C_PTR)",
    "pfft_timer" => "type(C_PTR)",
    "pfftf_timer" => "type(C_PTR)",
    "pfftl_timer" => "type(C_PTR)",
    "pfftq_timer" => "type(C_PTR)",
    "pfft_gctimer" => "type(C_PTR)",
    "pfftf_gctimer" => "type(C_PTR)",
    "pfftl_gctimer" => "type(C_PTR)",
    "pfftq_gctimer" => "type(C_PTR)",
    "void *" => "type(C_PTR)",
    "char *" => "type(C_PTR)",
    "double *" => "type(C_PTR)",
    "float *" => "type(C_PTR)",
    "long double *" => "type(C_PTR)",
    "float128__ *" => "type(C_PTR)",
    "pfft_complex *" => "type(C_PTR)",
    "pfftf_complex *" => "type(C_PTR)",
    "pfftl_complex *" => "type(C_PTR)",
    "pfftq_complex *" => "type(C_PTR)",
    "ptrdiff_t *" => "type(C_PTR)",
    );

# C->Fortran map of supported argument types
%arg_types = (
    "int" => "integer(C_INT), value",
    "unsigned" => "integer(C_INT), value",
    "size_t" => "integer(C_SIZE_T), value",
    "ptrdiff_t" => "integer(C_INTPTR_T), value",

    "pfft_r2r_kind" => "integer(C_FFTW_R2R_KIND), value",
    "pfftf_r2r_kind" => "integer(C_FFTW_R2R_KIND), value",
    "pfftl_r2r_kind" => "integer(C_FFTW_R2R_KIND), value",
    "pfftq_r2r_kind" => "integer(C_FFTW_R2R_KIND), value",

    "double" => "real(C_DOUBLE), value",
    "float" => "real(C_FLOAT), value",
    "long double" => "real(C_LONG_DOUBLE), value",
    "__float128" => "real(16), value",

    "pfft_complex" => "complex(C_DOUBLE_COMPLEX), value",
    "pfftf_complex" => "complex(C_DOUBLE_COMPLEX), value",
    "pfftl_complex" => "complex(C_LONG_DOUBLE), value",
    "pfftq_complex" => "complex(16), value",

    "pfft_plan" => "type(C_PTR), value",
    "pfftf_plan" => "type(C_PTR), value",
    "pfftl_plan" => "type(C_PTR), value",
    "pfftq_plan" => "type(C_PTR), value",
    "const pfft_plan" => "type(C_PTR), value",
    "const pfftf_plan" => "type(C_PTR), value",
    "const pfftl_plan" => "type(C_PTR), value",
    "const pfftq_plan" => "type(C_PTR), value",

    "pfft_gcplan" => "type(C_PTR), value",
    "pfftf_gcplan" => "type(C_PTR), value",
    "pfftl_gcplan" => "type(C_PTR), value",
    "pfftq_gcplan" => "type(C_PTR), value",
    "const pfft_gcplan" => "type(C_PTR), value",
    "const pfftf_gcplan" => "type(C_PTR), value",
    "const pfftl_gcplan" => "type(C_PTR), value",
    "const pfftq_gcplan" => "type(C_PTR), value",

    "pfft_timer" => "type(C_PTR), value",
    "pfftf_timer" => "type(C_PTR), value",
    "pfftl_timer" => "type(C_PTR), value",
    "pfftq_timer" => "type(C_PTR), value",
    "const pfft_timer" => "type(C_PTR), value",
    "const pfftf_timer" => "type(C_PTR), value",
    "const pfftl_timer" => "type(C_PTR), value",
    "const pfftq_timer" => "type(C_PTR), value",

    "pfft_gctimer" => "type(C_PTR), value",
    "pfftf_gctimer" => "type(C_PTR), value",
    "pfftl_gctimer" => "type(C_PTR), value",
    "pfftq_gctimer" => "type(C_PTR), value",
    "const pfft_gctimer" => "type(C_PTR), value",
    "const pfftf_gctimer" => "type(C_PTR), value",
    "const pfftl_gctimer" => "type(C_PTR), value",
    "const pfftq_gctimer" => "type(C_PTR), value",

    "const int *" => "integer(C_INT), dimension(*), intent(in)",
    "ptrdiff_t *" => "integer(C_INTPTR_T), dimension(*), intent(out)",
    "const ptrdiff_t *" => "integer(C_INTPTR_T), dimension(*), intent(in)",

    "const pfft_r2r_kind *" => "integer(C_FFTW_R2R_KIND), dimension(*), intent(in)",
    "const pfftf_r2r_kind *" => "integer(C_FFTW_R2R_KIND), dimension(*), intent(in)",
    "const pfftl_r2r_kind *" => "integer(C_FFTW_R2R_KIND), dimension(*), intent(in)",
    "const pfftq_r2r_kind *" => "integer(C_FFTW_R2R_KIND), dimension(*), intent(in)",

    "double *" => "real(C_DOUBLE), dimension(*), intent(out)",
    "const double *" => "real(C_DOUBLE), dimension(*), intent(in)",
    "float *" => "real(C_FLOAT), dimension(*), intent(out)",
    "const float *" => "real(C_FLOAT), dimension(*), intent(in)",
    "long double *" => "real(C_LONG_DOUBLE), dimension(*), intent(out)",
    "const long double *" => "real(C_LONG_DOUBLE), dimension(*), intent(in)",
    "__float128 *" => "real(16), dimension(*), intent(out)",
    "const __float128 *" => "real(16), dimension(*), intent(in)",

    "pfft_complex *" => "complex(C_DOUBLE_COMPLEX), dimension(*), intent(out)",
    "pfftf_complex *" => "complex(C_FLOAT_COMPLEX), dimension(*), intent(out)",
    "pfftl_complex *" => "complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out)",
    "pfftq_complex *" => "complex(16), dimension(*), intent(out)",

    "const pfft_complex *" => "complex(C_DOUBLE_COMPLEX), dimension(*), intent(in)",
    "const pfftf_complex *" => "complex(C_FLOAT_COMPLEX), dimension(*), intent(in)",
    "const pfftl_complex *" => "complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(in)",
    "const pfftq_complex *" => "complex(16), dimension(*), intent(in)",

    "void *" => "type(C_PTR), value",
    "FILE *" => "type(C_PTR), value",

    "const char *" => "character(C_CHAR), dimension(*), intent(in)",

    # Although the MPI standard defines this type as simply "integer",
    # if we use integer without a 'C_' kind in a bind(C) interface then
    # gfortran complains.  Instead, since MPI also requires the C type
    # MPI_Fint to match Fortran integers, we use the size of this type
    # (extracted by configure and substituted by the Makefile).
    "MPI_Comm" => "integer(\@C_MPI_FINT\@), value",
    "MPI_Comm *" => "integer(\@C_MPI_FINT\@), intent(out)"
    );

while (<>) {
    next if /^ *$/;
    if (/^ *extern +([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *\((.*)\) *$/) {
	$ret = &canonicalize_type($1);
	$name = $2;

	$args = $3;
	$args =~ s/^ *void *$//;

	$bad = ($ret ne "void") && !exists($return_types{$ret});	
	foreach $arg (split(/ *, */, $args)) {
	    $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
	    $argtype = &canonicalize_type($1);
	    $bad = 1 if !exists($arg_types{$argtype});
	}
	if ($bad) {
	    print "! Unable to generate Fortran interface for $name\n";
	    next;
	}

	# any function taking an MPI_Comm arg needs a C wrapper (grr).
	if ($args =~ /MPI_Comm/) {
	    $cname = $name . "_f03";
	}
	else {
	    $cname = $name;
	}

        # Since Fortran isn't case sensitive we must distiguish between get_N and get_n
        $name =~ s/get_n/get_Nos/;

	# Fortran has a 132-character line-length limit by default (grr)
	$len = 0;

	print "    "; $len = $len + length("    ");
	if ($ret eq "void") {
	    $kind = "subroutine"
	}
	else {
	    print "$return_types{$ret} ";
	    $len = $len + length("$return_types{$ret} ");
	    $kind = "function"
	}
	print "$kind $name("; $len = $len + length("$kind $name(");
	$len0 = $len;
	
	$argnames = $args;
	$argnames =~ s/([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) */$2/g;
	$comma = "";
	foreach $argname (split(/ *, */, $argnames)) {
            $argname = &canonicalize_name($argname);
	    if ($len + length("$comma$argname") + 3 > 132) {
		printf ", &\n%*s", $len0, "";
		$len = $len0;
		$comma = "";
	    }
	    print "$comma$argname";
	    $len = $len + length("$comma$argname");
	    $comma = ",";
	}
	print ") "; $len = $len + 2;

	if ($len + length("bind(C, name='$cname')") > 132) {
	    printf "&\n%*s", $len0 - length("$name("), "";
	}
	print "bind(C, name='$cname')\n";

	print "      import\n";
	foreach $arg (split(/ *, */, $args)) {
	    $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
	    $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
	    $ftype = $arg_types{$argtype};

	    # Various special cases for argument types:
	    if ($name =~ /_flops$/ && $argtype eq "double *") {
		$ftype = "real(C_DOUBLE), intent(out)" 
	    }
	    if ($name =~ /_execute/ && ($argname eq "ri" ||
					$argname eq "ii" || 
					$argname eq "in")) {
		$ftype =~ s/intent\(out\)/intent(inout)/;
	    }

	    print "      $ftype :: $argname\n"
	}

	print "    end $kind $name\n";
	print "    \n";
    }
}
