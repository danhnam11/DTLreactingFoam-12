// calculate CPU time for single steps
rhoTime = Time2 - Time1; 
UTime = Time3 - Time2; 
YEcalPropertiesTime = Time4 - Time3;
pTime = Time5 - Time4; 
updatePropertiesTime = Time8 - Time7; 
totalTime =Time8 - Time0;
total_preTime =Time8 - Time_pre;

// write the calculated CPU time to logFile 
logSummaryFile 
    << runTime.value() << ","
    << rhoTime << ","
    << UTime << ","
    << YEcalPropertiesTime << ","
    << pTime << ","
    << updatePropertiesTime << ","
    << totalTime << endl;


// show CPU time on the screen
/*
Info 
<< "Summary CPU time:" << '\n'
    << "Time" << tab
    << "rhoTime(s)" << tab
    << "UTime(s)" << tab
    << "YEcalPropertiesTime(s)" << tab
    << "pTime(s)" << tab
    << "updatePropertiesTime(s)" << tab
    << "totalTime(s)" << tab
    << "total_preTime(s)" << endl;

Info
    << runTime.value() << tab
    << rhoTime << tab
    << UTime << tab
    << YEcalPropertiesTime << tab
    << pTime << tab
    << updatePropertiesTime << tab
    << totalTime << tab
    << total_preTime << endl;
*/

// show CPU time measurement on the screen
Info << "Summary CPU time:" << endl;
Info << tab << "Time    = " << runTime.value() << endl;
Info << tab << "Step2   = " << UTime << " (s)," 
    << tab << 100*UTime/total_preTime << "[%], momemtum" << endl;
Info << tab << "Step345 = " << YEcalPropertiesTime << " (s)," 
    << tab << 100*YEcalPropertiesTime/total_preTime << "[%], Chem+Energy+Cal.Properties" << endl;
Info << tab << "Step678 = " << pTime << " (s)," 
    << tab << 100*pTime/total_preTime << "[%], Pressure" << endl;
Info << tab << "Step10  = " << updatePropertiesTime << " (s)," 
    << tab << 100*updatePropertiesTime/total_preTime << "[%], updateProperties" << endl;    
Info << tab << "total   = " << total_preTime << " (s)," 
    << tab << 100*total_preTime/total_preTime << "[%], totalTime" << endl;