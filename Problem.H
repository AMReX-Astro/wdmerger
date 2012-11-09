/* problem-specific Castro:: declarations go here */

void center_of_mass();

Real locWgtSumOneSide (const std::string& name, Real time, int idir, int side, int bdir);
 
Real volWgtSumOneSide (const std::string& name, Real time, int side, int bdir);

Real volProductSum (const std::string& name1, const std::string& name2, Real time);

Real locSquaredSum (const std::string& name, Real time, int idir);
