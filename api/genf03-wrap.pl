#!/usr/bin/perl -w
# Generate Fortran 2003 wrappers (which translate MPI_Comm from f2c) from
# function declarations of the form (one per line):
#     extern <type> pfft_<name>(...args...)
#     extern <type> pfft_<name>(...args...)
#     ...
# with no line breaks within a given function.  (It's too much work to
# write a general parser, since we just have to handle FFTW's header files.)
# Each declaration has at least one MPI_Comm argument.

sub canonicalize_type {
    my($type);
    ($type) = @_;
    $type =~ s/ +/ /g;
    $type =~ s/^ //;
    $type =~ s/ $//;
    $type =~ s/([^\* ])\*/$1 \*/g;
    $type =~ s/long double/R/;
    $type =~ s/ptrdiff_t/INT/;
    $type =~ s/pfftl_complex/C/;
    $type =~ s/pfftl_([A-Za-z0-9_]+)/PX($1)/;
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


while (<>) {
    next if /^ *$/;
    if (/^ *extern +([a-zA-Z_0-9 ]+[ \*]) *pfftl_([a-zA-Z_0-9]+) *\((.*)\) *$/) {
	$ret = &canonicalize_type($1);
	$name = $2;

	$args = $3;
	
	print "\n$ret PX(${name}_f03)(";

        $ret_comm = 0;
	$comma = "";
	foreach $arg (split(/ *, */, $args)) {
            $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
            $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
	    print $comma;
	    if ($argtype eq "MPI_Comm") {
		print "MPI_Fint f_$argname";
	    }
	    elsif ($argtype eq "MPI_Comm *") {
		print "MPI_Fint * f_$argname";
                $ret_comm = 1;
	    }
	    else {
		print "$argtype $argname";
	    }
	    $comma = ", ";
        }
	print ")\n{\n";

	print "  MPI_Comm ";
	$comma = "";
	foreach $arg (split(/ *, */, $args)) {
            $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
            $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
	    if (($argtype eq "MPI_Comm") || ($argtype eq "MPI_Comm *")) {
		print "$comma$argname";
		$comma = ", ";
	    }
        }
        print ";\n\n";

	foreach $arg (split(/ *, */, $args)) {
            $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
            $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
            if ($argtype eq "MPI_Comm") {
                print "  $argname = MPI_Comm_f2c(f_$argname);\n";
            }
        }

        $comma="";
        $argnames="";
	foreach $arg (split(/ *, */, $args)) {
            $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
            $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
            if ($argtype eq "MPI_Comm *") {
              $argnames = "$argnames$comma&$argname";
            }
            else {
              $argnames = "$argnames$comma$argname";
            }
            $comma=", ";
        }
# 	$argnames = $args;
# 	$argnames =~ s/([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) */$2/g;
	print "  ";
	print "$ret ret = " if ($ret ne "void");
	print "PX($name)($argnames);\n";
        if ($ret_comm) {
          foreach $arg (split(/ *, */, $args)) {
            $arg =~ /^([a-zA-Z_0-9 ]+[ \*]) *([a-zA-Z_0-9]+) *$/;
            $argtype = &canonicalize_type($1);
            $argname = &canonicalize_name($2);
            if ($argtype eq "MPI_Comm *") {
              print "  *f_$argname = MPI_Comm_c2f($argname);\n";
            }
          }
        }
	print "  return ret;\n" if ($ret ne "void");
	print "}\n";
    }
}
