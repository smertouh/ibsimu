#include <fstream>
#include <iomanip>
#include <limits>

#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include <string>
#include <sys/stat.h>
using namespace std;
#include <filesystem>



//const std::string bfieldfn = "/home/oleg/Desktop/ibsimu/simulations/2021/2021-01-22/SCHARGECOMPENSATION/By.dat";






class ForcedPot : public CallbackFunctorB_V {

public:

    ForcedPot() {}
    ~ForcedPot() {}

    virtual bool operator()( const Vec3D &x ) const {
        return( x[2] < 0.2e-3 && x[0]*x[0]+x[1]*x[1] > 6e-3*6e-3 );
    }
};


void simu( int argc, char **argv, double Vpuller, double Vdump, string current_working_dir )
{
    string sfile=argv[3];
    double comp=1.0e-3;
    double start = -3.0e-3;
    double h = 1e-3; //1e-3
    double sizereq[3] = { 140.0e-3,
                          140.0e-3, 
                          300.0e-3-start };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( -70.0e-3, -70.0e-3, start );
    Geometry geom( MODE_3D, meshsize, origo, h );

    
    double angle =  0.5*M_PI;
    std::cout << "angle = " << angle << "\n";
    Transformation T;
    T.translate( Vec3D( 0, 0, 0 ) );
    T.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    T.rotate_x( 0.5*M_PI );
    T.translate( Vec3D( 0, 0, h/50.0 ) );
    T.rotate_z( angle );
	string sbase=current_working_dir;//"/home/oleg/Desktop/ibsimu/simulations/2022";
    STLFile *fplasma = new STLFile(  sbase+"/PG.stl" );
    STLFile *fpuller = new STLFile( sbase+"/EG1.stl" );
    STLFile *fpuller2 = new STLFile( sbase+"/EG2.stl" );
    STLFile *fedump = new STLFile( sbase+"/AG.stl" );


    STLSolid *plasma = new STLSolid;
    plasma->set_transformation( T );
    plasma->add_stl_file( fplasma );
    geom.set_solid( 7, plasma );

    
    STLSolid *puller = new STLSolid;
    puller->set_transformation( T );
    puller->add_stl_file( fpuller );
    puller->add_stl_file( fpuller2 );
    geom.set_solid( 8, puller );

    STLSolid *edump = new STLSolid;
    edump->set_transformation( T );
    edump->add_stl_file( fedump );

    geom.set_solid( 9, edump );




    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vdump) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vpuller) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vdump) );



    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_Z, 0.2e-3 );
    solver.set_nsimp_initial_plasma( &initp );
    ForcedPot force;
    solver.set_forced_potential_volume( 0.0, &force );
    solver.set_gnewton( true );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    // Define magnetic field
    bool fout[3] = {true, true, true};
    string bfieldfn = sbase+"/By.dat";//+sfile;
    MeshVectorField bfield( MODE_3D, fout, 1.0e-3, 1.0, bfieldfn );
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO };
    bfield.set_extrapolation( bfldextrpl );

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 10000 ); //1000
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    // Suppress effects of magnetic field at volume where phi<1000. This is needed
    // because of erroneous field values close to the plasma electrode magnetic
    // steel (bug in Radia-3D).
    NPlasmaBfieldSuppression psup( epot, 1000.0 );
    pdb.set_bfield_suppression( &psup );

    double rho_h, rho_tot;
    size_t imax=15;
    for( size_t i = 0; i < imax; i++ ) { //15(1)
	
		ibsimu.message(1) << "Major cycle " << i << "/"<<imax<<"\n";
		ibsimu.message(1) << "-----------------------\n";
		
		
		if( i == 1 ) {
	    	std::vector<double> Ei, rhoi;
	    	Ei.push_back( 2.0 );
	    	rhoi.push_back( 0.5*rho_h );
	    	double rhop = rho_tot - rho_h*0.5;
            	solver.set_nsimp_plasma( rhop, 10.0, rhoi, Ei );
        	}
		
		solver.solve( epot, scharge_ave );
		efield.recalculate();
	
		double J = atof(argv[4]);//300.0;//22
        	pdb.clear(); 
		pdb.add_cylindrical_beam_with_energy( 50000, J, -1.0, 1.0, 
					      	2.0, 0.0, 1.0, 
					      	Vec3D(0,0,start),
					      	Vec3D(1,0,0), 
					      	Vec3D(0,1,0),
					      	18e-3 ); //50000(1000)
		rho_h = pdb.get_rhosum();
		pdb.add_cylindrical_beam_with_energy( 50000, J/6*10.0, -1.0, 1.0/1836.15, 
					      	2.0, 0.0, 1.0, 
					      	Vec3D(0,0,start),
					      	Vec3D(1,0,0), 
					      	Vec3D(0,1,0),
					      	18e-3 );//50000(1000)
        	pdb.iterate_trajectories( scharge, efield, bfield );
		rho_tot = pdb.get_rhosum();
	
        	
    		// Space charge compensation
        	for( uint32_t k = 0; k < scharge.size(2); k++ ) {       
            	double z = scharge.origo(2)+k*scharge.h();
            	for( uint32_t j = 0; j < scharge.size(1); j++ ) {   
                	//double y = scharge.origo(1)+j*scharge.h();
                	for( uint32_t i = 0; i < scharge.size(0); i++ ) {       
                    	//double x = scharge.origo(0)+i*scharge.h();
                    	if( z >= 90.0e-3 )
                        	scharge(i,j,k) *= comp;
                	}
            	}
        	}
    	// Space charge averaging
		if( i == 0 ) {
	    	scharge_ave = scharge;
		} else {
	    	double coef = 0.3;
	    	scharge *= coef;
	    	scharge_ave += scharge;
	    	scharge_ave *= (1.0/(1.0+coef));
		}
    	}
    //save yyp profile
    string outdir2 = sbase+sfile+"/";
    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_Z, geom.max(2)-geom.h(), 
                                   PARTICLE_DIAG_PLOT_SCATTER,
                                   DIAG_Y, DIAG_YP );
    pplotter2.set_size( 800, 600 );
    pplotter2.set_font_size( 16 );
    pplotter2.set_ranges( -0.03, -0.060, 0.03, 0.060 );
    pplotter2.set_emittance_ellipse( true );
    pplotter2.set_histogram_n( 151 );
    pplotter2.set_histogram_m( 151 );
    pplotter2.set_colormap_interpolation( INTERPOLATION_CLOSEST );
    pplotter2.plot_png( outdir2 + "emit_yyp.png" );
    pplotter2.export_data( outdir2 + "emit_yyp.txt" );
    //save xxp profile
    ParticleDiagPlotter pplotter1( geom, pdb, AXIS_Z, geom.max(2)-geom.h(), 
                                   PARTICLE_DIAG_PLOT_SCATTER,
                                   DIAG_X, DIAG_XP );
    pplotter1.set_size( 800, 600 );
    pplotter1.set_font_size( 16 );
    pplotter1.set_ranges( -0.02, -0.100, 0.05, 0.20 );
    pplotter1.set_emittance_ellipse( true );
    pplotter1.set_histogram_n( 151 );
    pplotter1.set_histogram_m( 151 );
    pplotter1.set_colormap_interpolation( INTERPOLATION_CLOSEST );
    pplotter1.plot_png( outdir2 + "emit_xxp.png" );
    pplotter1.export_data( outdir2 + "emit_xxp.txt" );
    
    //end profile
    geom.save( outdir2 +"geom.dat" );
    epot.save( outdir2 +"epot.dat" );
    pdb.save( outdir2 +"pdb.dat" );
    scharge.save(outdir2 +"scharge.dat");
    bfield.save(outdir2 +"bfield.dat");
}


int main( int argc, char **argv )
{
	char buff[FILENAME_MAX];
  	getcwd( buff, FILENAME_MAX );
  	string current_working_dir(buff);
	cout<<current_working_dir<<"\n";
	string sfile=argv[3];
	const char* chsfile;
	chsfile=sfile.c_str();
	mkdir(chsfile,0777);
    
    try {
	    string sfilelog=argv[3];
	ibsimu.set_message_output( sfilelog+"/ibsimu.txt" );//argv[3]+
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	double Vpuller = atof(argv[2]);
	double Vdump = atof(argv[1]);
	

	simu( argc, argv ,Vpuller, Vdump,current_working_dir);
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
