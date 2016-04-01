/*
A KBase module: MEModeling
*/

module MEModeling {
    /* 
        The workspace ID for a ContigSet data object.
        @id ws KBaseGenomes.Genome
    */
    typedef string ws_genome_id;

    /*
        A string representing a workspace name.
    */
    typedef string workspace_name;

    typedef structure {
	string report_name;
	string report_ref;
    } MEModelingResult;

    /*
        Build a Metabolic/Expression model
    */
    funcdef build_me_model(workspace_name workspace, ws_genome_id genome_id) returns (MEModelingResult) authentication required;
};