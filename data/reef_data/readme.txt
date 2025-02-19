Data was exported from the REEF database on September 15, 2023.

This dataset was provided to Brice Semmens. This dataset includes all data collected during REEF's Explore Baja Field Survey Trip in 2022.

**REEF requests that these data are not distributed to other parties without prior notification to REEF's Executive Director, and that these data are not put in to any data repositories such as GBIF, OBIS or IOOS, and that the geographic location data are not shared or published.**

The data should be cited as:
REEF. 2023. Reef Environmental Education Foundation. World Wide Web electronic publication. www.REEF.org, date of download (15 September 2023).

Files sent are: ExploreBaja22.txt, TEP-species.txt, and Bajageog.txt.

ExploreBaja22.txt file contain all fish sightings from Type 1 surveys conducted on the REEF Explore Baja Trip in 2022. The file contains 14,137 records, plus one row of column headings.

Each row is a fish sighting. Sightings from a given survey can be grouped using the form#.  Fields included are: 

1-	form #, unique # identifying a particular survey, multiple species sightings can be grouped together using this #),
2-	member number (confidential, not to be shared),
3-	experience of surveyor (Novice or Expert),
4-	geographic zone code (see the Bajageog.txt),
5-	site name,
6-	survey date,
7-	surface temperature in F.,
8-	bottom temperature in F.,
9-	bottom time, 
10-	start time,
11-	visibility,
12-	average depth (1 is snorkel, 2 is <10 feet, 3 is 10'-19', etc. up to 14 for 120'-129'),
13-	maximum depth (1 is snorkel, 2 is <10 feet, 3 is 10'-19', etc. up to 14 for 120'-129') - this information was added to the surveys starting in 2016,
14-	current (1-none, 2-week, 3-strong), and
15-	habitat type (there are 7 types of general habitats),
16-	REEF species code,
17-	REEF family code,
18-	REEF abundance codes (1- Single, 2- Few (2-10), 3- Many (11-100), 4- Abundant (>100))


Habitat codes are 
1 - Rock Boulder/Shale Reef is where a rocky or shale-like substrate protrudes from the sea floor, and may be covered in low algae (sea palms, feather boa, etc), coral, sponge, or other living organisms. 
2 - Wall is a vertical dropoff of over 20 feet that faces open water.
3 - Pinnacle is a large rock or group of rocks with steep sides that rises from deep water toward the surface.
4 - Sandy Bottom/Mud/Silt/Mangrove is where little to no rocks occur; bottom is composed of sand, sometimes mixed with silt or mud, may include mangroves.
5 - Open Ocean is a deep water area away from substrate where the bottom is not visible.
6 - Artificial includes shipwrecks, platforms, dumped debris or other artificially created habitats (incl. jetties).
7 - Rocky, Boulder, Sand Slope is a sloping drop-off where the bottom is a mixture of large and small rocks (typical habitat in the Galapagos).


Bajageog.txt contains a listing of the REEF Geographic Zone Codes surveyed on the Explore Baja trip, with site names and lat/lon if known. 


TEPspecies.txt contains a listing of all fish species in the REEF TEP database, including common name, scientific name, and REEF species ID code (as referenced in the fish file). 


The typical metrics that are used in analyses of REEF data are presence/absence, Sighting Frequency, Density Score, and Abundance Score.  These are calculated as -
Percent sighting frequency (%SF) for each species is the percentage of all dives in which the species was recorded.

Density score (D) for each species is a weighted average index based on the frequency of observations in different abundance categories.  Density score is calculated as: D= ((nSx1)+(nFx2)+(nMx3)+(nAx4)) / (nS + nF  + nM + nA), where nS, nF, nM, and nA represented the number of times each abundance category (Single, Few, Many, Abundant) was assigned for a given species.  Values range from 1 to 4.

Abundance Score is an estimate of abundance that accounts for non-sightings and is calculated as D x %SF, values range from 0 to 4.  

More information on these metrics, standard filters used on the dataset, and examples of how REEF data have been used can be found by reviewing the PDF papers and reports posted online at http://www.REEF.org.