# cs124_pa_02

command line instructions:<br />

$ ./strassen 0 dimension inputfile -- initial problem (strassen multiplicaiton)<br />

$ ./strassen tr p tests -- count triangles for the graph where edge appears with probability p, repeat tests times (default version tests = 5, p = 0.01)<br />

$ ./strassen spl dimension -- finds optimal split for given dimension (in case no argument dimension is given, runs over n = 2,4,..., 1024)

$ ./strassen test -- compares results by Strassen and trivial multiplication for n = 256; halts in case results don't match
