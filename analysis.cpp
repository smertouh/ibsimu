#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
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



using namespace std;



void simu( int argc, char **argv )
{
    string sfile=argv[6];
    string sbase="/home/oleg/Desktop/ibsimu/simulations/2021/2021-01-22/SCHARGECOMPENSATION/";
    string sbase2="/home/oleg/Desktop/ibsimu/simulations/2021/2021-01-22/SCHARGECOMPENSATION/traj/";
	string outdir2 = sbase+sfile+"/";
    std::ifstream is_geom( outdir2 +argv[1] );
    if( !is_geom.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[1] + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();
    geom.build_surface();

    std::ifstream is_epot( outdir2 +argv[2] );
    if( !is_epot.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[2] + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    std::ifstream is_pdb( outdir2 +argv[3] );
    if( !is_pdb.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[3] + "\'" ) );
    ParticleDataBase3D pdb( is_pdb, geom );
    is_pdb.close();
	
	MeshScalarField *scharge=NULL;
	if( argc >= 6 ){
	std::ifstream is_SCfield( outdir2 +argv[5] );
	MeshScalarField *mesh_SCfield = new MeshScalarField(is_SCfield);
	scharge= mesh_SCfield;
	}
	
    VectorField *bfield = NULL;
    if( argc >= 5 ) {
	bool fout[3] = {true, true, true};
	MeshVectorField *mesh_bfield = new MeshVectorField( MODE_3D, fout, 1.0e-3, 1.0, outdir2 +argv[4] );
	field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
					 FIELD_ZERO, FIELD_ZERO, 
					 FIELD_ZERO, FIELD_ZERO };
	mesh_bfield->set_extrapolation( bfldextrpl );
	bfield = mesh_bfield;
    }
    
	
    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
    if( bfield )
	plotter.set_bfield( bfield );
    plotter.set_trajdens( &tdens );
    if ( scharge )
    plotter.set_scharge( scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
    
	
	//new
	GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 0 );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_font_size( 16 );
    geomplotter.set_trajdens( &tdens );
    geomplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    geomplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );
	std::vector<double> eqlines;
	eqlines.push_back( 10 );
	eqlines.push_back( 1000 );
    eqlines.push_back( 2000 );
    eqlines.push_back( 4000 );
    eqlines.push_back( 6000 );
    eqlines.push_back( 6500 );
    eqlines.push_back( 7000.0 );
	eqlines.push_back( 8000.0 );
	eqlines.push_back( 10000.0 );
	eqlines.push_back( 20000.0 );
	eqlines.push_back( 21000.0 );
	eqlines.push_back( 22000.0 );
	eqlines.push_back( 23000.0 );
	eqlines.push_back( 24000.0 );
	eqlines.push_back( 25000.0 );
	eqlines.push_back( 25200.0 );
	eqlines.push_back( 25400.0 );
	eqlines.push_back( 25700.0 );
	eqlines.push_back( 85000.0 );
    geomplotter.set_eqlines_manual( eqlines );

    geomplotter.set_view( VIEW_ZX, -1 );
    geomplotter.plot_png( outdir2 + "geom_zx.png" );
    geomplotter.plot_png( sbase2+"zx/"+sfile + "geom_zx.png" );
    geomplotter.set_view( VIEW_ZY, -1 );
    geomplotter.plot_png( outdir2 + "geom_zy.png" );
    geomplotter.plot_png( sbase2+"zy/"+sfile + "geom_zy.png" );
    geomplotter.set_view_si( VIEW_XY, 5.0e-3 );
    geomplotter.plot_png( outdir2 + "geom_xy.png" );
    geomplotter.plot_png( sbase2+"xy/"+sfile + "geom_xy.png" );
	
	TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_X );
        diagnostics.push_back( DIAG_XP );
        diagnostics.push_back( DIAG_Y );
        diagnostics.push_back( DIAG_YP );
        pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diagnostics );
        Emittance emit_xxp( tdata(0).data(), tdata(1).data() );
        Emittance emit_yyp( tdata(2).data(), tdata(3).data() );

        // Output
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit_xxp.alpha() << " "
             << emit_xxp.beta() << " "
			 << emit_xxp.gamma() << " "
             << emit_xxp.epsilon() << " "
			 << emit_yyp.alpha() << " "
             << emit_yyp.beta() << " "
			 << emit_yyp.gamma() << " "
             << emit_yyp.epsilon() << "\n";
        dout.close();
        
        //save yyp profile
    
    //string outdir2 = sbase+sfile+"/";
    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_Z, geom.max(2)-geom.h(), 
                                   PARTICLE_DIAG_PLOT_SCATTER,
                                   DIAG_Y, DIAG_YP );
    pplotter2.set_size( 800, 600 );
    pplotter2.set_font_size( 16 );
    pplotter2.set_ranges( -0.005, -0.10, 0.005, 0.10 );
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
}


int main( int argc, char **argv )
{
    if( argc <= 3 ) {
	cerr << "Usage: analysis geom epot pdb <bfield>\n";
	exit( 1 );
    }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
