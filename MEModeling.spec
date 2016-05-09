/*
A KBase module: MEModeling
*/

module MEModeling {
    typedef structure {
		string model_ref;
		string workspace;
		string output_id;
    } BuildMEModelParams;

    typedef structure {
	string report_name;
	string report_ref;
	string me_model_ref;
    } MEModelingResult;

    /*
        Build a Metabolic/Expression model
    */
    funcdef build_me_model(BuildMEModelParams) returns (MEModelingResult) authentication required;
};