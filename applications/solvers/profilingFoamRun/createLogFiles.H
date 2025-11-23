
// create log File to assess computational time of each time step
OFstream logSummaryFile
(
    runTime.path()/("logCPUTime_" + runTime.userTimeName() + ".csv")
);

logSummaryFile 
//<< "variables=" << '\n'
    << "Time" << ","
    << "rhoTime(s)" << ","
    << "UTime(s)" << ","
    << "YEcalPropertiesTime(s)" << ","
    << "pTime(s)" << ","
    << "updatePropertiesTime(s)" << ","
    << "totalTime(s)" << endl;
