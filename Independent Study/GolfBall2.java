/**
	Graph plotting code taken from RealPendulumTraectoryStdents.java and ZombiesStudentVersion.java
	Author: J.-F. Briere & S. Bhatnagar
	Current version written: March 2019
*/

import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.util.Scanner;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import java.util.Arrays;

public class GolfBall2
{
	public static final boolean plotxVsy = true; // Only store and plot the angle vs time graph when set to true

    public static double accx (double vx, double vy, double v)
	{
		double Cd = 0.5;        // Drag Coefficient for non-dimpled Golf Ball
	    double So = 0.00006;       // Magnus coefficient
        double mass = 0.04593;  // kg
	    double area = 0.00143;  // m^2, radius = 21.335mm
        double rho = 1.184;     // kg/m^3
		double omegaX = 50.0;    // angular speed in x-direction rad/s 
		double omegaZ = 130.0;   // angular speed in z-direction rad/s
		
		if (v > 14.0)
		{
			Cd = 7.0/v;
		}
		
		double aDragx = ((0.5 * Cd * rho * area * (vx*vx))/mass);  // Formula for the acceleration due to drag in x-direction
		
		/**The Magnus force is given by the cross product of the angular velocity and linear velocity. 
	       The angular velocity is in the x and z-directions while the velocity is only in the x and y direction.
		   Using the formula for cross product, the Magnus force in the x-direction is given by Mx = (So/m)(-omegaZ*Vy)
		*/
	
	    double aMagnusx = ((So/mass) * (omegaZ * vy)); // x-component of the acceleration due to Magnus (So*Vy)
	
	    /** The differntial equation we are solving for tells us that the x-component of the acceleration of the ball 
		    is equal to the sum of the x components of the acceleration due to Drag and Magnus. Both of which are in negative direction
		*/
		
		double accx = (((-1) * (aDragx)) - (aMagnusx)); 
		return accx;
	}
	
	public static double accy (double vy, double vx, double vz, double v)
	{
		double Cd = 0.5;        
	    double So = 0.00006;        
        double mass = 0.04593;   
	    double area = 0.001432;
        double rho = 1.184;     
	    double omegaX = 50.0;    
		double omegaZ = 130.0;  
		double g = 9.81;        // accelerartion due to gravity, m/s^2
		
		if (v > 14.0)
		{
			Cd = 7.0/v;
		}
		
		double aDragy = ((0.5 * Cd * rho * area * (vy*vy))/mass); // Formula for the acceleration due to drag in y-direction 
		
		// After applying cross product formula the equation for Magnus in y-direction is given by My = (So/m) * (omegaZ * Vx)
		
		double aMagnusy = ((So/mass) * ((omegaZ * vx) - (omegaX * vz))); // y-component of the acceleration due to Magnus (So*Vx)
	    
		/** The y-component of the acceleration of the ball is equal to the sum of the y-components of the acceleration 
		    due to Magnus and Drag, minus the acceleration due to gravity g.
		*/
		
		double accy = (((-1) * (aDragy)) + (aMagnusy) - g);
		return accy;
	}
	
	public static double accz (double vz, double vy, double v)
	{
		double Cd = 0.5;        
	    double So = 0.00006;        
        double mass = 0.04593;   
	    double area = 0.001432;
        double rho = 1.184;     
	    double omegaX = 50.0;   
		double omegaZ = 130.0;
		
		if (v > 14.0)
		{
			Cd = 7.0/v;
		}
		
		double aDragz = ((0.5 * Cd * rho * area * (vz * vz))/mass); //Formula for the accelerartion due to drag in z-direction
		
		// After applying cross product formula the equation for Magnus in z-direction is given by Mz = (So/m) * (omegaX * Vy)
		
		double aMagnusz = ((So/mass) * (omegaX * vy)); // y-component of the acceleration due to Magnus (So*Vx)
	    
		/** The y-component of the acceleration of the ball is equal to the sum of the z-components of the acceleration 
		    due to Magnus and Drag
		*/
		
		double accz = (((-1) * aDragz) + (aMagnusz));
		return accz;
 	}	
	
	public static void main (String[] args)
	{
        double dt = 1.0e-3;     // time step, s
        double tMax = 12.0;     // max time 
	
	    int imax = (int)(tMax/dt);                    //index of arrays
	    double[] v = new double[imax];                //array for velocity of ball
		double[] vx = new double[imax];               //array for x-component of velocity of the ball
		double[] vy = new double[imax];               //array for y-component of velocity of the ball     
		double[] vz = new double[imax];               //array for z-component of velocity of the ball
	    double[] x = new double[imax];                //array for the x-position of the golf ball
	    double[] y = new double[imax];                //array for the y-position of the golf ball
		double[] z = new double[imax];                //array for the z-position of the golf ball
	    double theta = 0.00174533;                    //initial angle of ball with horizontal, 0.1 degrees
		double[] xOpt = new double[imax];			  //array for x-position of golf ball with optimum angle
		double[] yOpt = new double[imax];             //array for y-position of golf ball with optimum angle
		double[] zOpt = new double[imax];			  //array for z-position of golf ball with optimum angle
	    double[] x45 = new double[imax];              //array for the x-position of the golf ball with initial angle of 45 degrees
        double[] y45 = new double[imax];              //array for the y-position of the golf ball with initial angle of 45 degrees
		double[] z45 = new double[imax];              //array for the z-position of the golf ball with initial angle of 45 degrees
		double[] thetaG = new double[900];            //array for the initial angle of the ball to display on graph later, 900 vales since loop stops when angle is PI/2
		double[] xFinal = new double[900];            //array for final position of ball in x-position to display on graph later
		
		v[0] = 70.0;                     // initial speed of ball
		vx[0] = v[0] * Math.cos(theta);  // initial speed of x-component
		vy[0] = v[0] * Math.sin(theta);  // initial speed of y-component
		vz[0] = (130.0*0.021335);        // initial speed of z-component, omegaZ*radius
		x[0] = 0.0;                      // initial x-position of the ball
		y[0] = 0.0;                      // initial y-position of the ball
		z[0] = 0.0;                      // initial z-position of the ball
		
	    vx[1] = ((vx[0]) + ((accx(vx[0],vy[0], v[0])) * dt));                  // Euler's method to solve for x-component of velocity
		vy[1] = ((vy[0]) + ((accy(vy[0],vx[0], vz[0], v[0])) * dt));           // Euler's method to solve for y-component of velocity
		vz[1] = ((vz[0]) + ((accz(vz[0], vy[0], v[0])) * dt));                 // Euler's method to solve for z-component of velocity		
		v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); // Pythagorean theorem to find magnitude of velocity
		
		x[1] = (x[0] + (vx[1] * dt));  // Euler to get one value for position in x-direction
		y[1] = (y[0] + (vy[1] * dt));  // Euler to get one value for positon in y-direction 
		z[1] = (z[0] + (vz[1] * dt)); // Euler to get one value for position in z-direction  
		
		double xf = x[1]; //xf is the current final position of the ball
		
		/** Now we will start the loop to get the values for the position of the ball during its motion with an initial
		    angle of 0.1 degrees
		*/
		
		for (int i=2; y[i-2] >= 0.0; i++)
		{
		    vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * dt));              //Euler's method still used for components of velocity
		    vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * dt));
			vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * dt));
			
		    v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
		    x[i] = ((2*x[i-1]) - (x[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * (dt*dt)));         //Verlet's is now used to determine position of ball]
		    y[i] = ((2*y[i-1]) - (y[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * (dt*dt)));
			z[i] = ((2*z[i-1]) - (z[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * (dt*dt)));
	        
			xf = x[i]; //The value of xf is reassigned until the loop ends and we will take that final value of x[i] as the final position
		}
		
		double xMax = xf;           //xMax will be the current max distance of the ball
		double thetaMax = theta; // thetaMax will be the current optimum angle of the ball
		
		/** Now we start the loop to augment the angle to determine the optimum angle for max distance. The loop will start at 
		    0.2 degrees and be augmented by 0.1 degrees each time until the value of theta is equal to Pi.
		*/
		
		theta = 0.00349066;    // 0.2 degrees
		int ctr = 1;           // counter to assign values to xFinal and thetaG array later   
		
		while(theta < (Math.PI/2.0))
		{
		    vx[0] = v[0] * Math.cos(theta); //The initial velocity will always be 70 m/s, and the initial z-velocity will stay the same
			vy[0] = v[0] * Math.sin(theta);
	        x[0] = 0.0;
		    y[0] = 0.0;
			z[0] = 0.0;
			
		    /** The code will now execute the same sequence to find initial values for the x and y components before starting a
		        loop again to determine the rest of the ball's motion using Euler for the first part and Verlet for the second
		    */
		
            vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0])) * dt));
		    vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0], v[0])) * dt));
			vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0])) * dt)); 
			
		    v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1]));
		
		    x[1] = (x[0] + (vx[1] * dt));
		    y[1] = (y[0] + (vy[1] * dt));
			z[1] = (z[0] + (vz[1] * dt));
		
	        xf = x[1];
		
		    for (int i=2; y[i-2] >= 0.0; i++)
		    {
			    vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * dt));
		        vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * dt));
				vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * dt));
			
		        v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
		        x[i] = ((2*x[i-1]) - (x[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * (dt*dt)));
			    y[i] = (y[i-1] + (vy[i] * dt)); // Verlet's won't work here
			    z[i] = (z[i-1] + (vz[i] * dt));
				System.out.println(x[i] + "\t" + z[i]);
			    xf = x[i];
			}
	       
		    /** The final positions of x and initial angle of theta will be added to the xFinal and thetaG arrays respetively. That way they can be 			
			    displayed on a graph later.
			*/
			
			xFinal[ctr] = xf ;
			thetaG[ctr] = Math.toDegrees(theta);
			
			/** Before restarting the angle loop again, the new value of xf is compared to the current value of xMax. If the value
		        of xf is greater than xMax, xMax is changed to the value of xf and thetaMax is changed to the current value
			    of theta[0] which is the intial angle. If xf is less than xMax, nothing changes
	        */
			
			if (xf > xMax)
		    {
			    xMax = xf;
			    thetaMax = theta;
		    }
			
			theta += 0.00174533;
			ctr++;
		}	
		
		/** After finding the value of thetaMax, we need to reassign the values of the array since they are still equal to the 
		    last iteration of the loop. The sequence will still be exactly the same as before with Euler to start off and a for
		    loop with Verlet
		*/
		
		vx[0] = v[0] * Math.cos(thetaMax);
		vy[0] = v[0] * Math.sin(thetaMax);
		xOpt[0] = 0.0;
		yOpt[0] = 0.0;
		zOpt[0] = 0.0;
		
        vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0])) * dt));
		vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0])) * dt));
		vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0])) * dt));
		
		v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1]));  
		
		xOpt[1] = (xOpt[0] + (vx[1] * dt));
		yOpt[1] = (yOpt[0] + (vy[1] * dt));
		zOpt[1] = (zOpt[0] + (vz[1] * dt));
		
		for (int i=2; yOpt[i-1] >= 0.0; i++)
		{
			vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * dt));
		    vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * dt));
			vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * dt));
			
		    v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
		    xOpt[i] = (2*(xOpt[i-1])) - (xOpt[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * (dt*dt));
			yOpt[i] = (2*(yOpt[i-1])) - (yOpt[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * (dt*dt));
			zOpt[i] = (2*(zOpt[i-1])) - (zOpt[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * (dt*dt));
			
			xMax = xOpt[i];
	    }
		
	/** This section of the code was supposed to create new arrays for the values of the position of the golf ball with an intial	
	    angle of 45 degrees. That way there are two sets of arrays. One for optimum angle and max distance and the other set for
	    an intial angle of 45 degrees. And it will follow the exact same sequence as well
	*/

		theta = 0.785398; // 45 degrees
		vx[0] = v[0] * Math.cos(theta);
		vy[0] = v[0] * Math.sin(theta);
		x45[0] = 0.0;
		y45[0] = 0.0;
		z45[0] = 0.0;
		
	    vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0])) * dt));
		vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0])) * dt));
		vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0])) * dt));
		
		v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); 
		
		x45[1] = (x45[0] + (vx[1] * dt));
		y45[1] = (y45[0] + (vy[1] * dt));
		z45[1] = (z45[0] + (vz[1] * dt));
		double xMax45 = x45[1];
		
		for (int i=2; y45[i-1] >= 0.0; i++)
		{
		    vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * dt));
		    vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * dt));
			vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * dt));
			
		    v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
		    x45[i] = ((2*x45[i-1]) - (x45[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1])) * (dt*dt)));
			y45[i] = ((2*y45[i-1]) - (y45[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1])) * (dt*dt)));
			z45[i] = ((2*z45[i-1]) - (z45[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1])) * (dt*dt)));
			xMax45 = x45[i];
	    }
	    
	    /** And this section is all about graphing which is mostly cited so I don't think I have much to say here.
	        Hopefully these comments helped in some way
        */
		
	   ////////////////////////////////////////////////////////////////
		//1. Opening a file to store the x vs y data
		///////////////////////////////////////////////////////////////

		// Open log file to write the angle vs time data
  	  	String filename = "GolfSimulationOutput.txt";
		PrintWriter outputFile = null;
  	  	try
  	  	{
  	  	  outputFile = new PrintWriter(new FileOutputStream("x_vs_y.txt",false));
  	  	}
  	  	catch(FileNotFoundException e)
  	  	{
  	  	  System.out.println("File error.  Program aborted.");
  	  	  System.exit(0);
  	  	}
		
		//////////////////////////////////////////////////////////////////////////////
		// 3. Print and plot angle as a function of time
		//////////////////////////////////////////////////////////////////////////////
		if (plotxVsy)
	  	{
			//Storing info in a file to use later if needed
			for(int i=0;i<(imax-2);i+=10) // Only print every 10 points
			{
				outputFile.println(x[i] + "	" + y[i]); 
				outputFile.println(x45[i] + " " + y45[i]);
				outputFile.println(x[i] + " " + z[i]);
			}

	   		plot(thetaMax, xOpt, yOpt);//Method below
			plot(theta, x45,y45);
			plot(thetaMax, xOpt, zOpt); 
			plot2(thetaG, xFinal);
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// 4. Find and then print period as a function of amplitude
		//////////////////////////////////////////////////////////////////////////////
		System.out.println("x\t\t\ty");
		System.out.println("(m)\t\t\t(m)");
		System.out.println(xMax + "\t0");
		System.out.println(xMax45 + "\t0");
		
		///////////////////////////////////////////////////////////////////////
		// 5. Closing file
		///////////////////////////////////////////////////////////////////////
		outputFile.close();
	
	    System.out.printf("The angle in degrees to launch the golf ball at for maximum distance is %.3f\n",Math.toDegrees(thetaMax));
	    System.out.printf("The maximum distance travelled by the ball in meters is %.3f\n",xMax);	
		System.out.printf("The maximum distance travelled by the ball in meters when the angle is 45 degrees is %.3f\n",xMax45);
		
    }		
    
	////////////////////////////////////////////////////////////////////////////////////////////
  	// 7. Live graphing of angular postion as a function of time for amplitude = printamplitude
  	///////////////////////////////////////////////////////////////////////////////////////////
  	public static void plot(double angle, double[] x, double[] y)
  	{
  	  	Plot2DPanel plot1 = new Plot2DPanel();

    	// define the legend position
    	plot1.addLegend("SOUTH");

    	 // Add a line plot to the PlotPanel
        plot1.addLinePlot("Optimum Angle", x, y);
        plot1.setAxisLabel(0,"x (m)");
    	plot1.setAxisLabel(1,"y (m)");
		plot1.getAxis(0).setLabelPosition(0.5, -0.1);
        BaseLabel title1 = new BaseLabel("x v.s. y with Max angle "
    	  + Math.toDegrees(angle) + " degrees", Color.BLACK, 0.5, 1.1);
        title1.setFont(new Font("Courier", Font.BOLD, 14));
        plot1.addPlotable(title1);
        
        JFrame frame1 = new JFrame("Ball Trajectory");
    	frame1.setSize(500, 500);
    	frame1.setContentPane(plot1);
    	frame1.setVisible(true);
		frame1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  	}
	  	public static void plot2(double[] x, double[] y)
  	{
  	  	Plot2DPanel plot2 = new Plot2DPanel();

    	// define the legend position
    	plot2.addLegend("SOUTH");

    	// add a line plot to the PlotPanel
    	plot2.addLinePlot("GolfBall", x, y);
    	plot2.setAxisLabel(0,"Angle(rad)");
    	plot2.setAxisLabel(1,"x distance(m)");
    	BaseLabel title2 = new BaseLabel("Angle v.s. Max Distance ", Color.BLACK, 0.5, 1.1);

    	title2.setFont(new Font("Courier", Font.BOLD, 14));
    	plot2.addPlotable(title2);

    	// put the PlotPanel in a JFrame like a JPanel
    	JFrame frame2 = new JFrame("Angle vs Max Distance");
    	frame2.setSize(500, 500);
    	frame2.setContentPane(plot2);
    	frame2.setVisible(true);
		frame2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  	}
}