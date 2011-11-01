/*--------------------------------------------------------------------*/
/* 
	MAIN BODY 
*/
/*--------------------------------------------------------------------*/

DataSet myData = ReadDataFile ("Alignment.phylip");

DataSetFilter myFilter = CreateFilter (myData,1);

/* Force to accept rooted trees */
ACCEPT_ROOTED_TREES = 1;
 
/* Reading tree file in newick format, only one line */
/* Is this the only way? */ 
fscanf ("Tree_100_174.000.newick","String", Treestr);

UseModel(USE_NO_MODEL);
Tree chronogram = Treestr;


/* HarvestFrequencies (Freqs, myFilter, 1, 1, 1); */

             /* A     C     G      T*/
Freqs = {{0.25}{0.25}{0.25}{0.25}};


/* Model values */
/* JC all values equal 1 with equal base freqs*/
            
global          AC := 1;
global          AT := 1;
global          AG := 1;
global	        CG := 1;
global          CT := 1;
global          GT := 1;

                /* A        C        G          T */                           
GTRRateMatrix = {{*,       AC*t,     AG*t,     AT*t}
                 {AC*t,    *,        CG*t,     CT*t}
                 {AG*t,    CG*t,     *,        GT*t}
                 {AT*t,    CT*t,     GT*t,       * }};


Model GTR = (GTRRateMatrix, Freqs);

global siteRate 	= 1;
Tree   siteTree 	= Treestr;

chronoBranches = BranchLength(chronogram,-1);
chronoNames = BranchName(chronogram,-1);
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



