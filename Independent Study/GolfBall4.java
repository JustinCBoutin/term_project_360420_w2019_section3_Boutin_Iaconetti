import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.util.Scanner;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import java.util.Arrays;

public class GolfBall4
{
	public static final boolean plotzVsx = true; // Only store and plot the angle vs time graph when set to true
	public static double Cd = 0.5;                // Drag Coefficient
	public static double So = 0.00006;            // Magnus Coefficient
	public static final double mass = 0.04593;    // mass in kg
	public static final double radius = 0.021335; // radius of golf ball in meters
	public static final double area = 0.00143;    // area of golf ball in m^2
	public static final double rho = 1.184;       // density of air in kg/m^3
	public static final double g = 9.81;          // accelerartion due to gravity in m/s^2
	public static final double dt = 1.0e-3;       // time step in seconds
    public static final double tMax = 20.0;       // total amount of time
    public static double wx = 1.0;                // initial guess for optimum angular speed in x-direction
    public static double wz = 1.0;	              // initial guess for optimum angular speed in z-direction
	
	public static double accx (double vx, double vy, double v, double omegaZ)
	{
		if (v > 14.0)
		{
			Cd = 7.0/v; // Change in drag coefficient due to dimples on golf ball
		}
		
		double aDragx = ((0.5 * Cd * rho * area * (vx*vx))/mass);  // Formula for the acceleration due to drag in x-direction
		
		/**The Magnus force is given by the cross product of the angular velocity and linear velocity. 
	       The angular velocity is in the x and z-directions while the velocity is only in the x and y direction.
		   Using the formula for cross product, the Magnus force in the x-direction is given by Mx = (So/m)(-omegaZ*Vy)
		*/
	
	    double aMagnusx = ((So/mass) * (omegaZ * vy)); 
	
	    /** The differntial equation we are solving for tells us that the x-component of the acceleration of the ball 
		    is equal to the sum of the x components of the acceleration due to Drag and Magnus. Both of which are in negative direction
		*/
		
		double accx = (((-1) * (aDragx)) - (aMagnusx)); 
		return accx;
	}
	
	public static double accy (double vy, double vx, double vz, double v, double omegaX, double omegaZ)
	{
		if (v > 14.0)
		{
			Cd = 7.0/v;
		}
		
		double aDragy = ((0.5 * Cd * rho * area * (vy*vy))/mass); // Formula for the acceleration due to drag in y-direction 
		
		// After applying cross product formula the equation for Magnus in y-direction is given by My = (So/m) * (omegaZ * Vx) - (omegaX * Vz)
		
		double aMagnusy = ((So/mass) * ((omegaZ * vx) - (omegaX * vz))); 
	    
		/** The y-component of the acceleration of the ball is equal to the sum of the y-components of the acceleration 
		    due to Magnus and Drag, minus the acceleration due to gravity g.
		*/
		
		double accy = (((-1) * (aDragy)) + (aMagnusy) - g);
		return accy;
	}
	
	public static double accz (double vz, double vy, double v, double omegaX)
	{
		if (v > 14.0)
		{
			Cd = 7.0/v;
		}
		
		double aDragz = ((0.5 * Cd * rho * area * (vz * vz))/mass); //Formula for the accelerartion due to drag in z-direction
		
		// After applying cross product formula the equation for Magnus in z-direction is given by Mz = (So/m) * (omegaX * Vy)
		
		double aMagnusz = ((So/mass) * (omegaX * vy)); 
	    
		/** The y-component of the acceleration of the ball is equal to the sum of the z-components of the acceleration 
		    due to Magnus and Drag
		*/
		
		double accz = (((-1) * aDragz) + (aMagnusz));
		return accz;
 	}

	public static void main (String[] args)
	{
		int imax = (int)(tMax/dt);                    //index of arrays
	    double[] v = new double[imax];                //array for velocity of ball
		double[] vx = new double[imax];               //array for x-component of velocity of the ball
		double[] vy = new double[imax];               //array for y-component of velocity of the ball     
		double[] vz = new double[imax];               //array for z-component of velocity of the ball
	    double[] x = new double[imax];                //array for the x-position of the golf ball
	    double[] y = new double[imax];                //array for the y-position of the golf ball
		double[] z = new double[imax];                //array for the z-position of the golf ball
		double[] xOpt = new double[imax];			  //array for x-position of golf ball with optimum angle
		double[] yOpt = new double[imax];             //array for y-position of golf ball with optimum angle
		double[] zOpt = new double[imax];			  //array for z-position of golf ball with optimum angle
		
		double theta = 0.00174533;                    // initial angle of first golf shot
		double omegaX = 1.0;                          // initial angular speed in x-direction which will change over time
		double omegaZ = 1.0;                          // initial angular speed in z-direction which will change over time
		double xf = 0.0;                              // variable to hold final position of ball at the end of each loop in x-direction     
		double zf = 0.0;                              // variable to hold final position of ball at the end of each loop in z-direction
		double xMax = 300.0;                          // desired final position in x-direction which can be changed
		double zMax = 50.0;                           // desired final position in z-direction which can be changed
		
		v[0] = 70.0;                     // initial speed of ball
		vx[0] = v[0] * Math.cos(theta);  // initial speed of x-component
		vy[0] = v[0] * Math.sin(theta);  // initial speed of y-component
		vz[0] = (omegaZ*radius);        // initial speed of z-component, omegaZ*radius
		x[0] = 0.0;                      // initial x-position of the ball
		y[0] = 0.0;                      // initial y-position of the ball
		z[0] = 0.0;                      // initial z-position of the ball
		
	    vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0],omegaZ)) * dt));              // Euler's method to solve for x-component of velocity
		vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0],wx,omegaZ)) * dt));        // Euler's method to solve for y-component of velocity
		vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0],wx)) * dt));              // Euler's method to solve for z-component of velocity		
		v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); 			// Pythagorean theorem to find magnitude of velocity
		
		x[1] = (x[0] + (vx[1] * dt));  // Euler to get one value for position in x-direction
		y[1] = (y[0] + (vy[1] * dt));  // Euler to get one value for positon in y-direction 
		z[1] = (z[0] + (vz[1] * dt)); // Euler to get one value for position in z-direction  
		
		
		for (int i=2; y[i-1] > 0.0; i++)
		{
		    vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * dt));              //Euler's method still used for components of velocity
		    vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],wx,omegaZ)) * dt));
			vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1],wx)) * dt));
			
		    v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
		    x[i] = ((2*x[i-1]) - (x[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * (dt*dt)));         //Verlet's is now used to determine position of ball]
		    y[i] = ((2*y[i-1]) - (y[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],wx,omegaZ)) * (dt*dt)));
			z[i] = ((2*z[i-1]) - (z[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1],wx)) * (dt*dt)));
		
			//System.out.println(x[i] + "\t" + z[i]);
			xf = x[i];
			zf = z[i];
		}
		
		theta = 0.00349066; // initial angle incremented by 0.1 degrees 
		
		/** This next loop will increase the angular speed in the z-direction as the accelerartion in the x-direction is dependent 
		    on that componenet due to the magnus coefficient. The angular speed in the x-direction will be kept constant at the 
			initial guess of wx.
		**/
		
		while (xf < xMax)
		{ 
			vx[0] = v[0] * Math.cos(theta); 
			vy[0] = v[0] * Math.sin(theta); 
			vz[0] = (omegaZ*radius);        // initial speed of z-component, omegaZ*radius
			x[0] = 0.0;                      // initial x-position of the ball
			y[0] = 0.0;                      // initial y-position of the ball
			z[0] = 0.0;                      // initial z-position of the ball
		
		    vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0],omegaZ)) * dt));              // Euler's method to solve for x-component of velocity
			vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0],wx,omegaZ)) * dt));        // Euler's method to solve for y-component of velocity
			vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0],wx)) * dt));              // Euler's method to solve for z-component of velocity		
			v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); 			// Pythagorean theorem to find magnitude of velocity
		
			x[1] = (x[0] + (vx[1] * dt));  // Euler to get one value for position in x-direction
			y[1] = (y[0] + (vy[1] * dt));  // Euler to get one value for positon in y-direction 
			z[1] = (z[0] + (vz[1] * dt)); // Euler to get one value for position in z-direction  
		
		
			for (int i=2; y[i-1] > 0.0; i++)
			{
				vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * dt));              //Euler's method still used for components of velocity
				vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],wx,omegaZ)) * dt));
				vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1],wx)) * dt));
			
				v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
				x[i] = ((2*x[i-1]) - (x[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * (dt*dt)));         //Verlet's is now used to determine position of ball]
				y[i] = ((2*y[i-1]) - (y[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],wx,omegaZ)) * (dt*dt)));
				z[i] = ((2*z[i-1]) - (z[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1],wx)) * (dt*dt)));
					
				//System.out.println(x[i] + "\t" + z[i]);
				xf = x[i];
				zf = z[i];
			}
	
			theta += 0.00174533; // increase angle until desired final position in x is achieved
            omegaZ++;		     // increase omegaZ until desired position is reached
	    }
	
	    /** Now that there is a guess for the optimum omegaZ, the next loop will work to find the optimum omegaX. The omegaZ found will be used during
		    the loop, but is still subject to change as there is still a relation between omegaX and omegaZ. The omegaX will be incremented at the end
			of the loop each time, while the omegaZ must satisfy certain factors to be changed
		**/
	 
		while (zf < zMax)
		{
			vx[0] = v[0] * Math.cos(theta); 
			vy[0] = v[0] * Math.sin(theta); 
			vz[0] = (omegaZ*radius);        // initial speed of z-component, omegaZ*radius
			x[0] = 0.0;                      // initial x-position of the ball
			y[0] = 0.0;                      // initial y-position of the ball
			z[0] = 0.0;                      // initial z-position of the ball
		
			vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0],omegaZ)) * dt));              // Euler's method to solve for x-component of velocity
			vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0],omegaX,omegaZ)) * dt));        // Euler's method to solve for y-component of velocity
			vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0],omegaX)) * dt));              // Euler's method to solve for z-component of velocity		
			v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); 			// Pythagorean theorem to find magnitude of velocity
		
			x[1] = (x[0] + (vx[1] * dt));  // Euler to get one value for position in x-direction
			y[1] = (y[0] + (vy[1] * dt));  // Euler to get one value for positon in y-direction 
			z[1] = (z[0] + (vz[1] * dt)); // Euler to get one value for position in z-direction  
		
		
	     	for (int i=2; y[i-1] > 0.0; i++)
			{
				vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * dt));              //Euler's method still used for components of velocity
				vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],omegaX,omegaZ)) * dt));
				vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1],omegaX)) * dt));
			
				v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
				x[i] = ((2*x[i-1]) - (x[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * (dt*dt)));         //Verlet's is now used to determine position of ball]
				y[i] = ((2*y[i-1]) - (y[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],omegaX,omegaZ)) * (dt*dt)));
				z[i] = ((2*z[i-1]) - (z[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1],omegaX)) * (dt*dt)));
					
				//System.out.println(x[i] + "\t" + z[i]);
				xf = x[i];
				zf = z[i];
			}
			
			// if the final x-position is greater than the desired final position, then the omegaZ is decreased and the angle is increased
			if (xf > xMax)
			{
				omegaZ--;
				theta+=0.0017453;
			}
			
			// if the final x-position is less than the desired final position, then the omegaZ is increased and the angle is increased
			if (xf < xMax)
			{
				omegaZ++;
				theta+=0.0017453;
			}
			
			omegaX++;
		}
		
		// This last section fills the Optimum arrays to plot all the desired graphs at the end of the code, using the optimum spin in both directions
		
		xOpt[0] = 0.0;
		yOpt[0] = 0.0;
		zOpt[0] = 0.0;
		
		vx[1] = ((vx[0]) + ((accx(vx[0],vy[0],v[0],omegaZ)) * dt));              // Euler's method to solve for x-component of velocity
		vy[1] = ((vy[0]) + ((accy(vy[0],vx[0],vz[0],v[0],omegaX,omegaZ)) * dt));        // Euler's method to solve for y-component of velocity
		vz[1] = ((vz[0]) + ((accz(vz[0],vy[0],v[0],omegaX)) * dt));              // Euler's method to solve for z-component of velocity		
		v[1] = Math.sqrt((vx[1] * vx[1]) + (vy[1] * vy[1]) + (vz[1] * vz[1])); 			// Pythagorean theorem to find magnitude of velocity
		
		xOpt[1] = (xOpt[0] + (vx[1] * dt));  // Euler to get one value for position in x-direction
		yOpt[1] = (yOpt[0] + (vy[1] * dt));  // Euler to get one value for positon in y-direction 
		zOpt[1] = (zOpt[0] + (vz[1] * dt)); // Euler to get one value for position in z-direction  
		
		for (int i=2; yOpt[i-1] > 0.0; i++)
		{
			vx[i] = ((vx[i-1]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * dt));              //Euler's method still used for components of velocity
			vy[i] = ((vy[i-1]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],omegaX,omegaZ)) * dt));
			vz[i] = ((vz[i-1]) + ((accz(vz[i-1],vy[i-1],v[i-1],omegaX)) * dt));
		
			v[i] = Math.sqrt((Math.pow(vx[i],2.0)) + (Math.pow(vy[i],2.0)) + (Math.pow(vz[i],2.0)));
		
			xOpt[i] = ((2*xOpt[i-1]) - (xOpt[i-2]) + ((accx(vx[i-1],vy[i-1],v[i-1],omegaZ)) * (dt*dt)));         //Verlet's is now used to determine position of ball]
			yOpt[i] = ((2*yOpt[i-1]) - (yOpt[i-2]) + ((accy(vy[i-1],vx[i-1],vz[i-1],v[i-1],omegaX,omegaZ)) * (dt*dt)));
			zOpt[i] = ((2*zOpt[i-1]) - (zOpt[i-2]) + ((accz(vz[i-1],vy[i-1],v[i-1],omegaX)) * (dt*dt)));
					
			//System.out.println(x[i] + "\t" + z[i]);
			xf = xOpt[i];
			zf = zOpt[i];
		}
			
		////////////////////////////////////////////////////////////////
		//1. Opening a file to store the z vs x data
		///////////////////////////////////////////////////////////////

		// Open log file to write the angle vs time data
  	  	String filename = "GolfSimulationOutput.txt";
		PrintWriter outputFile = null;
  	  	try
  	  	{
  	  	  outputFile = new PrintWriter(new FileOutputStream("z_vs_x.txt",false));
  	  	}
  	  	catch(FileNotFoundException e)
  	  	{
  	  	  System.out.println("File error.  Program aborted.");
  	  	  System.exit(0);
  	  	}
		
		//////////////////////////////////////////////////////////////////////////////
		// 3. Print and plot angle as a function of time
		//////////////////////////////////////////////////////////////////////////////
		if (plotzVsx)
	  	{
			//Storing info in a file to use later if needed
			for(int i=0;i<(imax-2);i+=10) // Only print every 10 points
			{
				outputFile.println(xOpt[i] + "	" + zOpt[i]); 
			}

			plot(theta, xOpt, yOpt);
			Arrays.sort(xOpt);
			Arrays.sort(zOpt);
			plot(theta, zOpt, xOpt);//Method below
			
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// 4. Find and then print period as a function of amplitude
		//////////////////////////////////////////////////////////////////////////////
		System.out.println("x\t\t\tz");
		System.out.println("(m)\t\t\t(m)");
		System.out.println(xf + "\t" + zf);
		
		///////////////////////////////////////////////////////////////////////
		// 5. Closing file
		///////////////////////////////////////////////////////////////////////
		outputFile.close();
	
	    System.out.printf("The angle in degrees to launch the golf ball at for the needed distance is %.3f\n",Math.toDegrees(theta));
	    System.out.printf("The initial speed of the ball is %.3f\n",v[0]);	
		System.out.printf("The required spin on the ball in the x-direction is %.3f\n",omegaX);
		System.out.printf("The required spin on the ball in the z-direction is %.3f\n",omegaZ);
    }		
    
	////////////////////////////////////////////////////////////////////////////////////////////
  	// 7. Live graphing of angular postion as a function of time for amplitude = printamplitude
  	///////////////////////////////////////////////////////////////////////////////////////////
  	public static void plot(double angle, double[] z, double[] x)
  	{
  	  	Plot2DPanel plot1 = new Plot2DPanel();

    	// define the legend position
    	plot1.addLegend("SOUTH");

    	 // Add a line plot to the PlotPanel
        plot1.addLinePlot("Optimum Angle", z, x);
        plot1.setAxisLabel(0,"z (m)");
    	plot1.setAxisLabel(1,"x (m)");
		plot1.getAxis(0).setLabelPosition(0.5, -0.1);
        BaseLabel title1 = new BaseLabel("z v.s. x with angle "
    	  + Math.toDegrees(angle) + " degrees", Color.BLACK, 0.5, 1.1);
        title1.setFont(new Font("Courier", Font.BOLD, 14));
        plot1.addPlotable(title1);
        
        JFrame frame1 = new JFrame("Ball Trajectory");
    	frame1.setSize(500, 500);
    	frame1.setContentPane(plot1);
    	frame1.setVisible(true);
		frame1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  	}
}