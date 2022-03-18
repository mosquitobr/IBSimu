/* Mode 3D simulation of an time-of-fligth in a mass espectrometer model
 *
 * A E0 (eV), J=I/(r0*r0*M_PI) (A/m²), q/M beam propagates through a TOFL20.
 *
 * This is the first design to simulate 3D - TOFL20
 */


//access the the IBSimu library
#include <fstream>
#include "iomanip"
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "gtkwindow.hpp"
#include "iostream"
#include "sstream"


using namespace std;


// Declaration of global parameters
double Vdeteletrons = 37e2; // electrons detector
double Vgradepositiva = 5e2; // positive grid
double Vgradenegativa = -5e2;  // negative grid
double Vlente = -95e1; // lens
double Vdrift = -29.24e2; // drift tube
double Vdetions = -48e2; // ions detector
double f =  0.5;
double sc = 1/f;

double E01 = 0.012; // particles 1 starting energy (eV)
double E02 = 0.012; // particles 2 starting energy (eV)
double E03 = 0.012; // particles 3 starting energy (eV)

double Npart1 = 5e4; // particles 1 number quantities
double Npart2 = 5e4; // particles 2 number quantities
double Npart3 = 5e4; // particles 3 number quantities

double q = 1.0; // particle charge (in electron charges)
double m1 = 136.2; // particle mass1 (in atomic units)
double m2 = m1+1; // particle mass2 (in atomic units)
double m3 = m2+1; // particle mass3 (in atomic units)
double Tp = 1e-3; // praralel temperature (eV)
double Tt = 1e-4; // transverse temperature (eV)

double vq1 = 1.3884e4*sqrt(E01/m1); // particles 1 starting velocity (m/s)
double vq2 = 1.3884e4*sqrt(E02/m2); // particles 2 starting velocity (m/s)
double vq3 = 1.3884e4*sqrt(E02/m3); // particles 3 starting velocity (m/s)

double vp1 = 1.3884e4*sqrt(Tp/m1); // particles 1 starting paralel velocity (m/s)
double vp2 = 1.3884e4*sqrt(Tp/m2); // particles 2 starting paralel velocity (m/s)
double vp3 = 1.3884e4*sqrt(Tp/m3); // particles 3 starting paralel velocity (m/s)

double vt1 = 1.3884e4*sqrt(Tt/m1); // particles 1 starting transverse velocity (m/s)
double vt2 = 1.3884e4*sqrt(Tt/m2); // particles 2 starting transverse velocity (m/s)
double vt3 = 1.3884e4*sqrt(Tt/m3); // particles 3 starting transverse velocity (m/s)

double r0 = 5e-4; // ion beam radius (m)
double h = 1e-3; // thicness of cross overlap with particle radius beam X radius ionization source (m)

double I1 = (Npart1*vq1*1.602e-19)/h; // ion beam current 1
double I2 = (Npart2*vq2*1.602e-19)/h; // ion beam current 2
double I3 = (Npart2*vq3*1.602e-19)/h; // ion beam current 3

double J1 = I1/(r0*r0*M_PI); // beam linear current elect. density (A/m2)
double J2 = I2/(r0*r0*M_PI); // beam linear current elect. density (A/m2)
double J3 = I3/(r0*r0*M_PI); // beam linear current elect. density (A/m2)

// Simulation function
void simu( int argc, char **argv )
{
    // Create simulation box
    Geometry geom( MODE_3D, Int3D(sc*70+1,sc*70+1,sc*360+1), Vec3D(-0.035,-0.035,0), f*1e-3 );

    // Define solids by DXF file CAD draw
    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 2 );
    dxffile->read( "tofl203d.dxf" ); //file name DXF
    DXFSolid *s1 = new DXFSolid( dxffile, "deteletrons" );
    s1->scale( 1e-3 );
    s1->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "gradepositiva" );
    s2->scale( 1e-3 );
    s2->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "gradenegativa" );
    s3->scale( 1e-3 );
    s3->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "lente" );
    s4->scale( 1e-3 );
    s4->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 10, s4 );
    DXFSolid *s5 = new DXFSolid( dxffile, "drift" );
    s5->scale( 1e-3 );
    s5->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 11, s5 );
    DXFSolid *s6 = new DXFSolid( dxffile, "detions" );
    s6->scale( 1e-3 );
    s6->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 12, s6 );
    
    // Set boundary conditions
    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) ); // xmin
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) ); // xmax
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) ); // ymin
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) ); // ymax
    
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, Vdeteletrons) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Vgradepositiva) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, Vgradenegativa) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, Vlente) );
    geom.set_boundary( 11, Bound(BOUND_DIRICHLET, Vdrift) );
    geom.set_boundary( 12, Bound(BOUND_DIRICHLET, Vdetions) );
   
    geom.build_mesh();
    geom.build_surface();	
    geom.save( "tofgeom.dat" );
    
    // Construct the required fields
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                      FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                      FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    EpotBiCGSTABSolver solver( geom, 1.0e-4, 1000000, 1.0e-4, 10, true );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_surface_collision( false );

    for( size_t i = 0; i < 5; i++ ) {
    solver.solve( epot, scharge );
    efield.recalculate();
    //defining input particle beams sources
    pdb.clear();
	pdb.add_cylindrical_beam_with_velocity( Npart1, J1, q, m1, 
				     vq1, vp1, vt1, // velocity
				     Vec3D(0,0,0.081), //center
				     Vec3D(1,0,0), //dir1
				     Vec3D(0,1,0), //dir2
				     r0 );

	pdb.add_cylindrical_beam_with_velocity( Npart2, J2, q, m2, 
				     vq2, vp2, vt2, // velocity
				     Vec3D(0,0,0.081), //center
				     Vec3D(1,0,0), //dir1
				     Vec3D(0,1,0), //dir2
				     r0 );

	pdb.add_cylindrical_beam_with_velocity( Npart3, J3, q, m3, 
				     vq3, vp3, vt3, // velocity
				     Vec3D(0,0,0.081), //center
				     Vec3D(1,0,0), //dir1
				     Vec3D(0,1,0), //dir2
				     r0 );
	
         // Define electron beam
        pdb.add_cylindrical_beam_with_velocity( 1e4, 3.8, -1, 1.0/1836, 
				     //1.3884e4*sqrt(30*1836), 1.3884e4*sqrt(Tp*1836), 1.3884e4*sqrt(Tt*1836), 
				     1.88e6, 1.88e4, 5.949e4,
				     Vec3D(0,0,0.081), //center
				     -Vec3D(1,0,0), //-dir1
				     Vec3D(0,1,0), //dir2
				     r0 );
					      
	pdb.iterate_trajectories( scharge, efield, bfield );
	pdb.save("pdb.dat");
	
	//saving trajectory diagnostics
	std::vector<trajectory_diagnostic_e> diagnostics;
	diagnostics.push_back( DIAG_T );
	diagnostics.push_back( DIAG_X );
	diagnostics.push_back( DIAG_VX );
	diagnostics.push_back( DIAG_Y );
	diagnostics.push_back( DIAG_VY );
	diagnostics.push_back( DIAG_Z );
	diagnostics.push_back( DIAG_VZ );
	diagnostics.push_back( DIAG_MASS );
	diagnostics.push_back( DIAG_QM );
	diagnostics.push_back( DIAG_CURR );
	TrajectoryDiagnosticData tof;
	pdb.trajectories_at_plane( tof, AXIS_Z, 0.3549, diagnostics );
	tof.export_data( "tof.txt" );

	const TrajectoryDiagnosticColumn &time = tof(0);


	Histogram1D tof_histo( 150, time.data());
	vector<double> time_data = tof_histo.get_data();
	ofstream of_histo("tof_histo.txt" );
	for( uint32_t i = 0; i < tof_histo.n(); i++ ) {
	of_histo << tof_histo.coord(i) << " "
        	<< tof_histo(i) << "\n";
	}

	
	}
	//Write output file containing all particles sources
	string fnpout = "tof_out.txt";
	ofstream fileOut( fnpout.c_str() );
	for( size_t k = 0; k < pdb.size(); k++ ) {
	//Output coordinates
	Particle3D &pp = pdb.particle( k );
	fileOut << setw(12) << pp.IQ() << " ";
	fileOut << setw(12) << pp.m() << " ";
	for( size_t j = 0; j < 7; j ++ )
        	fileOut << setw(12) << pp(j) << " ";
	fileOut << "\n";
	}
	fileOut.close();
	
	
	// Set parameters for plotting
	GeomPlotter geomplotter( geom );
	geomplotter.set_size( 1200, 1200 );
	geomplotter.set_epot( &epot );
	geomplotter.set_particle_database( &pdb );
	geomplotter.set_view( VIEW_ZY, sc*35 );
	geomplotter.plot_png( "tofplot_zy.jpg" );
	geomplotter.set_view( VIEW_ZX, sc*35 );
	geomplotter.plot_png( "tofplot_zx.jpg" );
	geomplotter.set_view( VIEW_XY, sc*354.9 );
	geomplotter.plot_png( "tofplot_xy.jpg" );

	if ( true ) {
		MeshScalarField tdens( geom );
		pdb.build_trajectory_density_field( tdens );
		GTKPlotter plotter( &argc, &argv );
		plotter.set_geometry( &geom );
		plotter.set_epot( &epot );
		plotter.set_bfield( &bfield );
		plotter.set_efield( &efield );
		plotter.set_scharge( &scharge_ave );
		plotter.set_trajdens( &tdens );
		plotter.set_particledatabase( &pdb );
		plotter.new_geometry_plot_window();
		plotter.run();
	}
}

int main( int argc, char **argv )
{
	try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 ) ;
	ibsimu.set_thread_count ( 4 );
		simu( argc, argv );
	} catch ( Error e ) {
		e.print_error_message( ibsimu.message( 0 ) );
		exit( 1 );
	}
    
	return ( 0 );

}
