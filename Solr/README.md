1. Installation

Follow the 'Solr Quick Start' link on the page:

http://lucene.apache.org/solr/resources.html#solr-version-control

to Download/Unpack/Launch Solr on your Solr server machine, if locally, it will be running at:

http://localhost:8983/solr/#/

2. After Solr is installed and up running, create the core directories for each core.
e.g., to create cores compounds and reactions, follow the steps below:

1) Under the solr installation where 'bin/solr start [stop]' can be fired:

-------------------------------------------------------------------------
mkdir my_solr_install_dir/server/solr/compounds
mkdir -p my_solr_install_dir/server/solr/compounds/conf
mkdir -p my_solr_install_dir/server/solr/compounds/data
mkdir my_solr_install_dir/server/solr/reactions
mkdir -p my_solr_install_dir/server/solr/reactions/conf
mkdir -p my_solr_install_dir/server/solr/reactions/data

-------------------------------------------------------------------------

2) Inside each of the above conf folders, create/copy-paste the solrconfig.xml, schema.xml 
and schema_types.xml. I copy-pasted these files from the corresponding dirs of the model_compound
and model_reaction of the patricbrc github repo (https://github.com/PATRIC3/patric_solr).

Note: From my email conversation with Hyunseung Yoo, 
"I should still define the type 'string_ci' inside schema_types.xml, except replace the 
'solr.TextDocValueField' with 'solr.SortableTextField'."

At this point the data folders are empty.

3) In the browser where you can see the Solr Admin UI, (i.e., http://localhost:8983/solr/#/), 
'Core Admin'-->'Add Core'. Set the instanceDir and dataDir by entering the absolute paths of the dirs
created in Step 1):
e.g., for compounds:
        instanceDir= my_solr_install_dir/server/solr/compounds
        dataDir= my_solr_install_dir/server/solr/compounds/data
      for reactions:
        instanceDir= my_solr_install_dir/server/solr/reactions
        dataDir= my_solr_install_dir/server/solr/reactions/data

3. Load data into Solr
At the shell go to where the bin/solr command is and run the following commands to load data
into the corresponding core we have just created in the browser UI. (Note: if you skipped step 3),
you will get "HTTP ERROR 404" (Not Found) error.

localhost:my_solr_install_dir qzhang$ bin/post -c core_name path-to-data/datafile.json

--------------When data loading returns error (it happens!)--------------

localhost:my_solr_install_dir qzhang$ bin/post -c reactions path-to-biochem_data/reactions.json

SimplePostTool version 5.0.0
Posting files to [base] url http://localhost:8983/solr/reactions/update...
Entering auto mode. File endings considered are xml,json,jsonl,csv,pdf,doc,docx,ppt,pptx,xls,xlsx,odt,odp,ods,ott,otp,ots,rtf,htm,html,txt,log
POSTing file reactions.json (application/json) to [base]/json/docs
SimplePostTool: WARNING: Solr returned an error #400 (Bad Request) for url: http://localhost:8983/solr/reactions/update/json/docs
SimplePostTool: WARNING: Response: {
  "responseHeader":{
    "status":400,
    "QTime":3424},
  "error":{
    "metadata":[
      "error-class","org.apache.solr.common.SolrException",
      "root-error-class","java.lang.NumberFormatException"],
    "msg":"ERROR: [doc=rxn13120] Error adding field 'deltag'='null' msg=For input string: \"null\"",
    "code":400}}
SimplePostTool: WARNING: IOException while reading response: java.io.IOException: Server returned HTTP response code: 400 for URL: http://localhost:8983/solr/reactions/update/json/docs
1 files indexed.
COMMITting Solr index changes to http://localhost:8983/solr/reactions/update...
Time spent: 0:00:06.774

-------------------------------------------------------------------------

The above error was caused by the .json file which contains data value of "null" as a string with quotes
instead of just null. So replacing the string "null" values with null without the quotes should solve the problem.

I also may need to replace /"rxn\d+" (space)/ (ignore the forward slashes in the regexpression) with nothing for the whole rections.json file in addition to "null" to null.
In the meantime, I renamed the field "source" in reactions.json to "data_source" because in the schema.xml file, the copyfield
is using "source" vs "dest".  To avoid confusion, I changed the "source" in reaction.json to "data_source".

Similarly, I need to replace /"cpd\d+" (space)/ (ignore the forward slashes in the regexpression) with nothing for the whole compounds.json file in addition to "null" to null.

And finally I should change the reactions/compounts.json files from a large object (i.e., {...}) into an array
of objects (i.e., [{...}, ..., {...}]) by replacing the first brace '{' with '[' and the last '}' with ']'.

-------------------------------------------------------------------------
localhost:my_solr_install_dir qzhang$ bin/post -c reactions path-to-biochem_data/reactions.json

java -classpath localhost:my_solr_install_dir/dist/solr-core-7.4.0.jar -Dauto=yes -Dc=reactions -Ddata=files org.apache.solr.util.SimplePostTool path-to-biochem_data/reactions.json
SimplePostTool version 5.0.0
Posting files to [base] url http://localhost:8983/solr/reactions/update...
Entering auto mode. File endings considered are xml,json,jsonl,csv,pdf,doc,docx,ppt,pptx,xls,xlsx,odt,odp,ods,ott,otp,ots,rtf,htm,html,txt,log
POSTing file Reactions.json (application/json) to [base]/json/docs
1 files indexed.
COMMITting Solr index changes to http://localhost:8983/solr/reactions/update...
Time spent: 0:00:06.997

Kept getting the "msg":"copyField dest :'data_source_str' is not an explicit field and doesn't match a dynamicField."
when tried to:
localhost:my_solr_install_dir qzhang$ bin/post path-to-biochem_data/reactions.json -c reactions
after the reactions.json file has been reformatted by replacing the /"rxn\d+" (space)/ with nothing and "null" with null.
Could not find a solution by editing the schema.xml file directly and even by deleting existing managed-schema file. 
Eventually I restored the schema.xml back to the previous one without the field of "source" or "data_source", and then
added the field via the Solr web API.  Then it worked.
Looking into the folder 'my_solr_install_dir/server/solr/reactions/conf', I found that the
newly added field 'data_source' and copyfield 'data_source' appear in the file managed-schema instead of in schema.xml.

After that, the following commands ran without a problem:

-------------------------------------------------------------------------

localhost:my_solr_install_dir qzhang$ bin/post path-to-biochem_data/reactions.json -c reactions

java -classpath localhost:my_solr_install_dir/dist/solr-core-7.4.0.jar -Dauto=yes -Dc=reactions -Ddata=files org.apache.solr.util.SimplePostTool path-to-biochem_data/reactions.json
SimplePostTool version 5.0.0
Posting files to [base] url http://localhost:8983/solr/reactions/update...
Entering auto mode. File endings considered are xml,json,jsonl,csv,pdf,doc,docx,ppt,pptx,xls,xlsx,odt,odp,ods,ott,otp,ots,rtf,htm,html,txt,log
POSTing file reactions.json (application/json) to [base]/json/docs
1 files indexed.
COMMITting Solr index changes to http://localhost:8983/solr/reactions/update...
Time spent: 0:00:09.381

-------------------------------------------------------------------------

At this point, refresh the browser Solr Admin UI, I see 34711 reaction records in the reactions core. 
and I can query the reactions indexes.

=========================================================================

4. The above procedures are for running Solr on my Mac laptop.  To run similar procedure on the megilen cloud,
where the cloud version of Solr (7) was installed, every step is the same, except the core creation (via the browser UI) requires the elevation.xml file in addition to the other .xml files.  So I simply copy-pasted the eleveation.xml file from
the megilen_solr_installation_dir/solr/server/solr/configsets/sample_techproducts_configs/conf/
directory into each of the compounds and reactions' conf directory.  Then run the following:

-------------------------------------------------------------------------

megilen-cloud-server:megilen_solr_installation_dir modelseedSolrUser$ bin/post -c compounds path-to-biochem_data/compounds.json

java -classpath megilen_solr_installation_dir/solr/dist/solr-core-7.4.0-SNAPSHOT.jar -Dauto=yes -Dc=compounds -Ddata=files org.apache.solr.util.SimplePostTool path-to-biochem_data/compounds.json
SimplePostTool version 5.0.0
Posting files to [base] url http://localhost:8983/solr/compounds/update...
Entering auto mode. File endings considered are xml,json,jsonl,csv,pdf,doc,docx,ppt,pptx,xls,xlsx,odt,odp,ods,ott,otp,ots,rtf,htm,html,txt,log
POSTing file compounds.json (application/json) to [base]/json/docs
1 files indexed.
COMMITting Solr index changes to http://localhost:8983/solr/compounds/update...
Time spent: 0:00:13.602

megilen-cloud-server:megilen_solr_installation_dir modelseedSolrUser$ bin/post -c reactions path-to-biochem_data/reactions.json
java -classpath megilen_solr_installation_dir/dist/solr-core-7.4.0-SNAPSHOT.jar -Dauto=yes -Dc=reactions -Ddata=files org.apache.solr.util.SimplePostTool path-to-biochem_data/reactions.json
SimplePostTool version 5.0.0
Posting files to [base] url http://localhost:8983/solr/reactions/update...
Entering auto mode. File endings considered are xml,json,jsonl,csv,pdf,doc,docx,ppt,pptx,xls,xlsx,odt,odp,ods,ott,otp,ots,rtf,htm,html,txt,log
POSTing file Reactions.json (application/json) to [base]/json/docs
1 files indexed.
COMMITting Solr index changes to http://localhost:8983/solr/reactions/update...
Time spent: 0:00:12.632

----------------------------------------------------------------------

5. To access the local Solr instance from local modelseed-ui
Inside the config.js of modelseed-ui, add an entry of 'local_solr_url: "http://0.0.0.0:8983/solr/"' to this.services.

Then in file 'modelseed-ui/app/services/biochem.js' and module 'Biochem', find the following line and set the solr_endpoint:

var solr_endpoint = config.services.local_solr_url;

After that, make sure wherever module 'Biochem' is injected and a call to 'Biochem.get*' is changed to'Biochem.get*_solr'.

Then you can see the biochemistry data changes in your local Solr reflected in modelseed-ui.
