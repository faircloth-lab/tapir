/*--------------------------------------------------------------------*/
/* 
	MAIN BODY 
*/
/*--------------------------------------------------------------------*/

//fscanf(stdin, "String", Alignment_File)
//fprintf(stdout, Alignment_File)
DataSet myData = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter myFilter = CreateFilter (myData,1);

/* Force to accept rooted trees */
ACCEPT_ROOTED_TREES = 1;
 
fscanf(stdin, "String", FILE_TREE);
//FILE_TREE = "Euteleost.tree";

    fscanf (FILE_TREE, REWIND, "Raw", treeString);

    treeStringPattern = treeString$"^#NEXUS";
    if (treeStringPattern[0] >= 0) {
        ExecuteCommands (treeString);
        if (IS_TREE_PRESENT_IN_DATA == 0) {
            fprintf (stdout, "\nThis NEXUS file doesn't contain a valid tree block");
            return 1;
        }
        // always take the first tree found
        treeString = NEXUS_FILE_TREE_MATRIX[0][1];
    } else {
        treeStringPattern = treeString$"\\(";
        if (treeStringPattern[0] < 0) {
            fprintf (stdout, "\nThis doesn't seem to be a valid Newick string file. Can't find the opening parenthesis.\nHad:", treeString).
            return 1;
        } else {
            parenCounter = 1;
            strlength = Abs (treeString);
            cp = treeStringPattern[0] + 1;
            while ( cp < strlength && parenCounter ) {
                cpc = treeString[cp];
                if (cpc == "(") {
                    parenCounter = parenCounter + 1;
                } else {
                    if (cpc == ")") {
                        parenCounter = parenCounter - 1;
                    }
                }
                cp = cp + 1;
            }

            if (parenCounter) {
                fprintf (stdout, "\nThis doesn't seem to be a valid Newick string file. Can't match the parentheses.\nHad:", treeString).
                return 1;
            }

            treeStringPattern = treeStringPattern[0];
            treeString = treeString[treeStringPattern][cp - 1];
        }
    }

Treestr = treeString;

UseModel(USE_NO_MODEL);
Tree chronogram = Treestr;


HarvestFrequencies (Freqs, myFilter, 1, 1, 1);

             /* A     C     G      T*/
// Freqs = {{0.25}{0.25}{0.25}{0.25}};
// fscanf(stdin, "Matrix", Freqs);

/* Model values */
/* JC all values equal 1 with equal base freqs*/

fscanf(stdin, "Matrix", Rates);
//Rates = {{0.9628}{0.5839}{1}{0.3660}{1.8756}{0.5105}};

global          AC := Rates[0];
global          AT := Rates[1];
global          AG := Rates[2];
global	        CG := Rates[3];
global          CT := Rates[4];
global          GT := Rates[5];

/*
global          AC := 1;
global          AT := 1;
global          AG := 1;
global	        CG := 1;
global          CT := 1;
global          GT := 1;
*/

                /* A        C        G          T */                           
GTRRateMatrix = {{*,       AC*t,     AG*t,     AT*t}
                 {AC*t,    *,        CG*t,     CT*t}
                 {AG*t,    CG*t,     *,        GT*t}
                 {AT*t,    CT*t,     GT*t,       * }};


Model GTR = (GTRRateMatrix, Freqs);

global siteRate 	= 1;
Tree   siteTree 	= Treestr;

chronoBranches = BranchLength(chronogram,-1);
chronoNames	   = BranchName (chronogram,-1);
chronoLength   = 0;
	
for (i=0; i< Columns(chronoBranches); i=i+1) 	{
	   chronoLength = chronoLength + chronoBranches[i];
	   ExecuteCommands ("siteTree." +  chronoNames[i] + ".t:=" + chronoBranches[i] + "*siteRate;");
}

fprintf (stdout, "\nChronogram length (time units): ", chronoLength, "\n");
fprintf (stdout, "\nBase freqs: ", Freqs, "\n");
fprintf (stdout, "\nSubstitution matrix values:\n",
                              "\t\tAC: ", AC, "\n",
                              "\t\tAG: ", AG, "\n",
                              "\t\tAT: ", AT, "\n",
                              "\t\tCG: ", CG, "\n",
                              "\t\tCT: ", CT, "\n",
                              "\t\tGT: ", GT, "\n\n\n",);

doneSites    = {myFilter.unique_sites,3}; /* Getting the different site patterns */
fullSites    = {myFilter.sites,3};
GetDataInfo    (dupInfo, myFilter); /*DupInfo is a matrix with the different categories/patterns of sites, duplicated sites info*/
alreadyDone  = {};

for (siteCount = 0; siteCount < myFilter.sites; siteCount = siteCount+1) {
			siteMap = dupInfo[siteCount];
			if (alreadyDone[siteMap] == 0) {
				filterString = "" + siteCount;
				DataSetFilter siteFilter = CreateFilter (myData,1,filterString);
				
                /* to avoid that 'siteRate' variable was stuck
                at the value estimated from the previous site */
                siteRate = 1;
                
                LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
                Optimize (site_res, siteLikelihood);
				alreadyDone[siteMap]  = 1;
        
                siteBranches = BranchLength (siteTree,-1);
                siteLength = 0;
	            for (i=0; i< Columns(siteBranches); i=i+1) 	{
            	   siteLength = siteLength + siteBranches[i];
                }
   				siteRate = siteLength/chronoLength;


				doneSites[siteMap][0] = siteLength;
                doneSites[siteMap][1] = siteRate;
				doneSites[siteMap][2] = site_res[1][0];/*Loglikelihood*/
                
			}
			ReportSite (siteCount, siteMap);				
}	



fprintf (stdout, "THE END\n");


/*------------------------------------------------------------------------*/

function ReportSite (siteI, siteM)
{
	fullSites[siteI][0] = doneSites[siteM][0]; 
	fullSites[siteI][1] = doneSites[siteM][1];
	fullSites[siteI][2] = doneSites[siteM][2];


	fprintf (stdout, " Site ", Format(siteI+1,4,0),
                     " Total subst = ", Format(fullSites[siteI][0],7,4), " subst,",
					 " Rate = ", Format(fullSites[siteI][1],7,4), " subst/time,",
					 " Log(L) ", Format(fullSites[siteI][2],7,4),"\n");		
					 
	return 0;
}

/*------------------------------------------------------------------------*/



