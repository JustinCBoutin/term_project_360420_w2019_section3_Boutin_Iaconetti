/**
	Graph plotting code taken from RealPendulumTraectoryStdents.java
	Author: J.-F. Briere
	Current version written: March 2019
*/

import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.util.Scanner;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import java.util.Arrays;

public class GolfBall
{
	public static final boolean plotxVsy = true; // Only store and plot the angle vs time graph when set to true

    public static double accx (double v, double theta)
	{
		double Cd = 0.5;        // Drag Coefficient for non-dimpled Golf Ball
	    double So = 0.22;       // (Magnus coefficient x omega)/mass 
        double mass = 0.04593;  // kg
	    double area = 0.00143;  // m^2, radius = 21.335mm
        double airDen = 1.184;     // density of air kg/m^3
		
        double aDrag = ((Cd * airDen * area * (v*v))/mass);        // Formula for the acceleration due to drag Cd*p*area*V^2
		double aDragx = (aDrag * Math.cos(theta));              // x-component of the acceleration due to drag 
		
		/**The Magnus force is given by the cross product of the angular velocity and linear velocity. 
	       If the angular velocity is purely in the z-direction, the Magnus force only acts in the x and y direction.
		   And after applying cross product formula the equation for Magnus in x-direction is given by Mx = So(-wz * Vy) 
		*/
	
	    double aMagnusx = (So * (v * Math.sin(theta))); // x-component of the acceleration due to Magnus (So*Vy)
	
	    /** The differntial equation we are solving for tells us that the x-component of the acceleration of the ball 
		    is equal to the sum of the x components of the acceleration due to Drag and Magnus.
		*/
		
		double accx = (((-1) * (aDragx)) - (aMagnusx)); 
		return accx;
	}
	
	public static double accy (double v, double theta)
	{
		double Cd = 0.5;        
	    double So = 0.212;        
        double mass = 0.0459;   
	    double area = 0.001432;
        double airDen = 1.168;     
	    double g = 9.81;        // accelerartion due to gravity, m/s^2
		
		double aDrag = ((Cd * airDen * area * (v*v))/mass);
		double aDragy = (aDrag * Math.sin(theta)); // y-component of the accelerartion due to drag
		
		// After applying cross product formula the equation for Magnus in y-direction is given by My = So(wz * Vx)
		
		double aMagnusy = (So * (v * Math.cos(theta))); // y-component of the acceleration due to Magnus (So*Vx)
	    
		/** The y-component of the acceleration of the ball is equal to the sum of the y-components of the acceleration 
		    due to Magnus and Drag, minus the acceleration due to gravity g.
		*/
		
		double accy = (((-1) * (aDragy)) + (aMagnusy) - g);
		return accy;
	}
	
	public static void main (String[] args)
	{
        double dt = 1.0e-3;     // time step, s
        double tMax = 12.0;     // max time 
	
	    int imax = (int)(tMax/dt);                    //index of arrays
	    double[] v = new double[imax];                //array for the speed of the golf ball
	    double[] x = new double[imax];                //array for the x-position of the golf ball
	    double[] y = new double[imax];                //array for the y-position of the golf ball
	    double[] theta = new double[imax];            //array for the angle of the ball with respect to the horizontal
	    double[] x45 = new double[imax];              //array for the x-position of the golf ball with initial angle of 45 degrees
        double[] y45 = new double[imax];              //array for the y-position of the golf ball with initial angle of 45 degrees
		double[] thetaG = new double[900];            //array for the initial angle of the ball to display on graph later, 900 vales since loop stops when angle is PI/2
		double[] xFinal = new double[900];            //array for final position of ball in x-position to display on graph later
		
		double[] xOpt = new double[imax];
		double[] yOpt = new double[imax];
		
		v[0] =70.0;           // initial speed of the ball  
		theta[0] = 0.00174533; // initial angle of the ball = 0.1 degrees
		x[0] = 0.0;            // initial x-position of the ball
		y[0] = 0.0;            // initial y-position of the ball
		
	    double vx = ((v[0] * Math.cos(theta[0])) + ((accx(v[0],theta[0])) * dt)); // Euler's method to solve for x-component of velocity
		double vy = ((v[0] * Math.sin(theta[0])) + ((accy(v[0],theta[0])) * dt)); // Euler's method to solve for y-component of velocity
		v[1] = Math.sqrt((vx * vx) + (vy * vy));                // Pythagorean theorem to find magnitude of velocity
		theta[1] = Math.atan(vy/vx);                                              // arctan of vy/vx to find the new angle of velocity vector with respect to horizontal
		
		x[1] = (x[0] + (vx * dt)); // Euler to get one value for position in x-direction for Verlet's later
		y[1] = (y[0] + (vy * dt)); // Euler to get one value for positon in y-direction for Verlet's later
		
		double xf = x[1]; //xf is the current final position of the ball
		
		/** Now we will start the loop to get the values for the position of the ball during its motion with an initial
		    angle of 0.1 degrees
		*/
		
		for (int i=2; y[i-1] > 0.0; i++)
		{
		    vx = ((v[i-1] * Math.cos(theta[i-1])) + ((accx(v[i-1],theta[i-1])) * dt)); //Euler's method still used for components of velocity
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + ((accy(v[i-1],theta[i-1])) * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    x[i] = ((x[i-1]) + (vx * dt)); //Verlet's is now used to determine position of ball
		    y[i] = ((y[i-1]) + (vy * dt));
			
	        xf = x[i+1]; //The value of xf is reassigned until the loop ends and we will take that final value of x[i] as the final position
		}
		
		double xMax = xf;           //xMax will be the current max distance of the ball
		double thetaMax = theta[0]; // thetaMax will be the current optimum angle of the ball
		
		/** Now we start the loop to augment the angle to determine the optimum angle for max distance. The loop will start at 
		    0.2 degrees and be augmented by 0.1 degrees each time until the value of theta is equal to Pi.
		*/
		
		theta[0] = 0.00349066; // 0.2 degrees
		int ctr = 1;           // counter to assign values to xFinal and thetaG array later   
		
		while(theta[0] < (Math.PI/2.0))
		{
		    v[0] = 70.0; //The initial velocity will always be 70 m/s
	        x[0] = 0.0;
		    y[0] = 0.0;
			
		    /** The code will now execute the same sequence to find initial values for the x and y components before starting a
		        loop again to determine the rest of the ball's motion using Euler for the first part and Verlet for the second
		    */
		
            vx = ((v[0] * Math.cos(theta[0])) + ((accx(v[0],theta[0])) * dt));
		    vy = ((v[0] * Math.sin(theta[0])) + ((accy(v[0],theta[0])) * dt));
		    v[1] = Math.sqrt((vx*vx) + (vy*vy));
		    theta[1] = Math.atan(vy/vx);
		
		    x[1] = (x[0] + (vx * dt));
		    y[1] = (y[0] + (vy * dt));
		
	        xf = x[1];
		
		    for (int i=2; y[i-1] > 0.0; i++)
		    {
			    vx = ((v[i-1] * Math.cos(theta[i-1])) + ((accx(v[i-1],theta[i-1])) * dt));
		        vy = ((v[i-1] * Math.sin(theta[i-1])) + ((accy(v[i-1],theta[i-1])) * dt));
		        v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		        theta[i] = Math.atan(vy/vx);
		
		        x[i] = (x[i-1] + (vx * dt));
			    y[i] = (y[i-1] + (vy * dt));
			
			    xf = x[i];
			}
	       
		    /** The final positions of x and initial angle of theta will be added to the xFinal and thetaG arrays respetively. That way they can be 			
			    displayed on a graph later.
			*/
			
			xFinal[ctr] = xf ;
			thetaG[ctr] = Math.toDegrees(theta[0]);
			//System.out.println(xFinal[ctr]);
			//System.out.println(thetaG[ctr]);
			
			/** Before restarting the angle loop again, the new value of xf is compared to the current value of xMax. If the value
		        of xf is greater than xMax, xMax is changed to the value of xf and thetaMax is changed to the current value
			    of theta[0] which is the intial angle. If xf is less than xMax, nothing changes
	        */
			
			if (xf > xMax)
		    {
			    xMax = xf;
			    thetaMax = theta[0];
		    }
			
			theta[0]+=0.00174533;
			ctr++;
		}	
		
		/** After finding the value of thetaMax, we need to reassign the values of the array since they are still equal to the 
		    last iteration of the loop. The sequence will still be exactly the same as before with Euler to start off and a for
		    loop with Verlet
		*/
		
		v[0] = 70.0;
		theta[0] = thetaMax;
		xOpt[0] = 0.0;
		yOpt[0] = 0.0;
		
        vx = ((v[0] * Math.cos(theta[0])) + ((accx(v[0],theta[0])) * dt));
		vy = ((v[0] * Math.sin(theta[0])) + ((accy(v[0],theta[0])) * dt));
		v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		theta[1] = Math.atan(vy/vx);
		
		xOpt[1] = (x[0] + (vx * dt));
		yOpt[1] = (y[0] + (vy * dt));
		
		for (int i=2; yOpt[i-1] > 0.0; i++)
		{
			vx = ((v[i-1] * Math.cos(theta[i-1])) + ((accx(v[i-1],theta[i-1])) * dt));
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + ((accy(v[i-1],theta[i-1])) * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    xOpt[i] = (2*(xOpt[i-1])) - (xOpt[i-2]) + ((accx(v[i-1],theta[i-1])) * (dt*dt));
			yOpt[i] = (2*(yOpt[i-1])) - (yOpt[i-2]) + ((accy(v[i-1],theta[i-1])) * (dt*dt));
			xMax = xOpt[i];	
	    }
		
	/** This section of the code was supposed to create new arrays for the values of the position of the golf ball with an intial	
	    angle of 45 degrees. That way there are two sets of arrays. One for optimum angle and max distance and the other set for
	    an intial angle of 45 degrees. And it will follow the exact same sequence as well
	*/

        v[0] = 70.0;
		theta[0] = 0.785398; // 45 degrees
		x45[0] = 0.0;
		y45[0] = 0.0;
		
	    vx = ((v[0] * Math.cos(theta[0])) + ((accx(v[0],theta[0])) * dt));
		vy = ((v[0] * Math.sin(theta[0])) + ((accy(v[0],theta[0])) * dt));
		v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		theta[1] = Math.atan(vy/vx);
		
		x45[1] = (x[0] + (vx * dt));
		y45[1] = (y[0] + (vy * dt));
		double xMax45 = x45[1];
		
		for (int i=2; y45[i-1] > 0.0; i++)
		{
		    vx = ((v[i-1] * Math.cos(theta[i-1])) + ((accx(v[i-1],theta[i-1])) * dt));
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + ((accy(v[i-1],theta[i-1])) * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    x45[i] = ((x45[i-1]) + (vx * dt));
			y45[i] = ((y45[i-1]) + (vy * dt));
			xMax45 = x45[i];
	    }
	    
	    /** And this section is all about graphing which is pretty much just copy paste so I don't think I have much to say here.
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
			}

	   		plot(thetaMax, xOpt, yOpt); //Method below
			plot(0.7853981634, x45, y45);
			plot2(thetaG, xFinal);
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// 4. Find and then print period as a function of amplitude
		//////////////////////////////////////////////////////////////////////////////
		System.out.println("x	        y");
		System.out.println("(m)		(m)");
		System.out.println(xMax + "  0 ");
		System.out.println(xMax45 + " 0 ");
		
		///////////////////////////////////////////////////////////////////////
		// 5. Closing file
		///////////////////////////////////////////////////////////////////////
		outputFile.close();
	
	    System.out.printf("The angle in degrees to launch the golf ball at for maximum distance is %.3f\n",Math.toDegrees(thetaMax));
	    System.out.printf("The maximum distance travelled by the ball in meters is %.3f\n",xMax);	
		System.out.printf("The maximum distance travelled by the ball in meters when the angle is 45 degrees is %.3f\n",xMax45);
		
    }		
    
	///////////////////////////////////////////////////////////////////////
  	// 7. Live graphing of angular postion as a function of time for amplitude = printamplitude
  	///////////////////////////////////////////////////////////////////////
  	public static void plot(double angle, double[] x, double[] y)
  	{
  	  	Plot2DPanel plot1 = new Plot2DPanel();

    	// define the legend position
    	plot1.addLegend("SOUTH");

    	// add a line plot to the PlotPanel
    	plot1.addLinePlot("GolfBall", x, y);
    	plot1.setAxisLabel(0,"x (m)");
    	plot1.setAxisLabel(1,"y (m)");
    	BaseLabel title1 = new BaseLabel("x v.s. y for max angle "
    	  + Math.toDegrees(angle) + " degrees", Color.BLACK, 0.5, 1.1);

    	title1.setFont(new Font("Courier", Font.BOLD, 14));
    	plot1.addPlotable(title1);

    	// put the PlotPanel in a JFrame like a JPanel
    	JFrame frame1 = new JFrame("Ball Trajectory");
    	frame1.setSize(800, 800);
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
    	frame2.setSize(800, 800);
    	frame2.setContentPane(plot2);
    	frame2.setVisible(true);
		frame2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  	}
}
