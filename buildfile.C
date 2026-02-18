{
    //making and setting a directory for libraries
    gSystem->mkdir("build", true);
    gSystem->SetBuildDir("build", true);

    //compile classes macros
    gSystem->CompileMacro("CYL.cxx", "k");
    gSystem->CompileMacro("Point.cxx", "k");
    gSystem->CompileMacro("VTX.cxx", "k");
    gSystem->CompileMacro("Particle.cxx", "k");

    //compile program macros
    gSystem->CompileMacro("simulation.C", "k");
    gSystem->CompileMacro("reconstruction.C", "k");

}
