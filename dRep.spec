/*
A KBase module: dRep
*/

module dRep {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef dereplicate(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
