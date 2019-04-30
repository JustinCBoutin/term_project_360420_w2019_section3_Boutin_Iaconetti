import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.util.Scanner;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import java.util.Arrays;

//Useless comment
public class GolfBall
{
	public static final boolean plotxVsy = true; // Only store and plot the angle vs time graph when set to true
	
	public static void main (String[] args)
	{
        double Cd = 0.1;     // Drag Coefficient for non-dimpled Golf Ball
	    double So = 0.25;       // (So*w)/mass
        double mass = 0.0459;   // kg
	    double area = 0.001432; // m^2
        double rho = 1.168;     // kg/m^3
	    double g = 9.8;         // m/s^2
	    double dt = 1.0e-3;     // time step, s
        double tMax = 15.0;     // max time
	    
		int imax = (int)(tMax/dt); //index of arrays
		double[] adrag = new double[imax];
		double[] amagnusx = new double[imax];
		double[] amagnusy = new double[imax];
		double[] v = new double[imax];
		double[] x = new double[imax];
		double[] y = new double[imax];
		double[] theta = new double[imax];
		
		v[0] = 70.0;
		theta[0] = 0.00174533;
		x[0] = 0.0;
		y[0] = 0.0;
		
        adrag[0] = ((Cd * rho * area * (v[0]*v[0]))/mass);
		double adragx = (adrag[0] * Math.cos(theta[0]));
		double adragy = (adrag[0] * Math.sin(theta[0]));
		  
		amagnusx[0] = (So * (v[0] * Math.sin(theta[0])));
		amagnusy[0] = (So * (v[0] * Math.cos(theta[0])));
		  
		double accx = (((-1) * (adragx)) - (amagnusx[0]));
		double accy = (((-1) * (adragy)) + (amagnusy[0]) - g);
		
		double vx = ((v[0] * Math.cos(theta[0])) + (accx * dt));
		double vy = ((v[0] * Math.sin(theta[0])) + (accy * dt));
		v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		theta[1] = Math.atan(vy/vx);
		
		x[1] = (x[0] + (vx * dt));
		y[1] = (y[0] + (vy * dt));
		
		double xf = x[1];
		
		for (int i=2; y[i-1] > 0.0; i++)
		{
			if(v[i-1] < 14.0)
			{
				Cd = 0.5;
			}
			
			else
			{
				Cd = 7.0/v[i-1];
			}
			
			adrag[i-1] = ((Cd * rho * area * (v[i-1]*v[i-1]))/mass);
			adragx = (adrag[i-1] * Math.cos(theta[i-1]));
		    adragy = (adrag[i-1] * Math.sin(theta[i-1]));
			
			amagnusx[i-1] = (So * (v[i-1] * Math.sin(theta[i-1])));
		    amagnusy[i-1] = (So * (v[i-1] * Math.cos(theta[i-1])));
			
			accx = (((-1) * (adragx)) - (amagnusx[i-1]));
		    accy = (((-1) * (adragy)) + (amagnusy[i-1]) - g);
			
			vx = ((v[i-1] * Math.cos(theta[i-1])) + (accx * dt));
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + (accy * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    x[i] = (2*(x[i-1])) - (x[i-2]) + (accx*(dt*dt));
			y[i] = (2*(y[i-1])) - (y[i-2]) + (accy*(dt*dt));
			
			xf = x[i];
		}
		
		double xMax = xf;
		double thetaMax = theta[0];
		
		for(theta[0] = 0.00349066; theta[0]<Math.PI; theta[0]+=0.00174533)
		{
		    v[0] = 70.0;
	        x[0] = 0.0;
		    y[0] = 0.0;
		    Cd = 7.0/v[0];
		
            adrag[0] = ((Cd * rho * area * (v[0]*v[0]))/mass);
		    adragx = (adrag[0] * Math.cos(theta[0]));
		    adragy = (adrag[0] * Math.sin(theta[0]));
		  
		    amagnusx[0] = (So * (v[0] * Math.sin(theta[0])));
		    amagnusy[0] = (So * (v[0] * Math.cos(theta[0])));
		  
		    accx = (((-1) * (adragx)) - (amagnusx[0]));
		    accy = (((-1) * (adragy)) + (amagnusy[0]) - g);
		
		    vx = ((v[0] * Math.cos(theta[0])) + (accx * dt));
		    vy = ((v[0] * Math.sin(theta[0])) + (accy * dt));
		    v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[1] = Math.atan(vy/vx);
		
		    x[1] = (x[0] + (vx * dt));
		    y[1] = (y[0] + (vy * dt));
		
	        xf = x[1];
		
		    for (int i=2; y[i-1] > 0.0; i++)
		    {
			    if(v[i-1] < 14.0)
			    {
				    Cd = 0.5;
			    }
			
			    else
			    {
				    Cd = 7.0/v[i-1];
			    }
			
			    adrag[i-1] = ((Cd * rho * area * (v[i-1]*v[i-1]))/mass);
			    adragx = (adrag[i-1] * Math.cos(theta[i-1]));
		        adragy = (adrag[i-1] * Math.sin(theta[i-1]));
			
			    amagnusx[i-1] = (So * (v[i-1] * Math.sin(theta[i-1])));
		        amagnusy[i-1] = (So * (v[i-1] * Math.cos(theta[i-1])));
			
			    accx = (((-1) * (adragx)) - (amagnusx[i-1]));
		        accy = (((-1) * (adragy)) + (amagnusy[i-1]) - g);
			
			    vx = ((v[i-1] * Math.cos(theta[i-1])) + (accx * dt));
		        vy = ((v[i-1] * Math.sin(theta[i-1])) + (accy * dt));
		        v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		        theta[i] = Math.atan(vy/vx);
		
		        x[i] = (2*(x[i-1])) - (x[i-2]) + (accx*(dt*dt));
		        y[i] = (y[i-1]) + (vy*dt);
			
			    xf = x[i];
		    }
		
		    if (xf > xMax)
		    {
			    xMax = xf;
			    thetaMax = theta[0];
		    }
	    }	
		
		v[0] = 70.0;
		theta[0] = thetaMax;
		x[0] = 0.0;
		y[0] = 0.0;
		Cd = 7.0/v[0];
		
        adrag[0] = ((Cd * rho * area * (v[0]*v[0]))/mass);
		adragx = (adrag[0] * Math.cos(theta[0]));
		adragy = (adrag[0] * Math.sin(theta[0]));
		  
		amagnusx[0] = (So * (v[0] * Math.sin(theta[0])));
		amagnusy[0] = (So * (v[0] * Math.cos(theta[0])));
		  
		accx = (((-1) * (adragx)) - (amagnusx[0]));
		accy = (((-1) * (adragy)) + (amagnusy[0]) - g);
		
		vx = ((v[0] * Math.cos(theta[0])) + (accx * dt));
		vy = ((v[0] * Math.sin(theta[0])) + (accy * dt));
		v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		theta[1] = Math.atan(vy/vx);
		
		x[1] = (x[0] + (vx * dt));
		y[1] = (y[0] + (vy * dt));
		
		for (int i=2; y[i-1] > 0.0; i++)
		{
			if(v[i-1] < 14.0)
			{
				Cd = 0.5;
			}
			
			else
			{
				Cd = 7.0/v[i-1];
			}
			
			adrag[i-1] = ((Cd * rho * area * (v[i-1]*v[i-1]))/mass);
			adragx = (adrag[i-1] * Math.cos(theta[i-1]));
		    adragy = (adrag[i-1] * Math.sin(theta[i-1]));
			
			amagnusx[i-1] = (So * (v[i-1] * Math.sin(theta[i-1])));
		    amagnusy[i-1] = (So * (v[i-1] * Math.cos(theta[i-1])));
			
			accx = (((-1) * (adragx)) - (amagnusx[i-1]));
		    accy = (((-1) * (adragy)) + (amagnusy[i-1]) - g);
			
			vx = ((v[i-1] * Math.cos(theta[i-1])) + (accx * dt));
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + (accy * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    x[i] = (2*(x[i-1])) - (x[i-2]) + (accx*(dt*dt));
			y[i] = (2*(y[i-1])) - (y[i-2]) + (accy*(dt*dt));
	    }
		
	/*	double[] x45 = new double[imax];
		double[] y45 = new double[imax];
		
		v[0] = 70.0;
		theta[0] = 0.785398;
		x45[0] = 0.0;
		y45[0] = 0.0;
		Cd = 7.0/v[0];
		
        adrag[0] = ((Cd * rho * area * (v[0]*v[0]))/mass);
		adragx = (adrag[0] * Math.cos(theta[0]));
		adragy = (adrag[0] * Math.sin(theta[0]));
		  
		amagnusx[0] = (So * (v[0] * Math.sin(theta[0])));
		amagnusy[0] = (So * (v[0] * Math.cos(theta[0])));
		  
		accx = (((-1) * (adragx)) - (amagnusx[0]));
		accy = (((-1) * (adragy)) + (amagnusy[0]) - g);
		
		vx = ((v[0] * Math.cos(theta[0])) + (accx * dt));
		vy = ((v[0] * Math.sin(theta[0])) + (accy * dt));
		v[1] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		theta[1] = Math.atan(vy/vx);
		
		x45[1] = (x[0] + (vx * dt));
		y45[1] = (y[0] + (vy * dt));
		
		for (int i=2; y[i-1] > 0.0; i++)
		{
			if(v[i-1] < 14.0)
			{
				Cd = 0.5;
			}
			
			else
			{
				Cd = 7.0/v[i-1];
			}
			
			adrag[i-1] = ((Cd * rho * area * (v[i-1]*v[i-1]))/mass);
			adragx = (adrag[i-1] * Math.cos(theta[i-1]));
		    adragy = (adrag[i-1] * Math.sin(theta[i-1]));
			
			amagnusx[i-1] = (So * (v[i-1] * Math.sin(theta[i-1])));
		    amagnusy[i-1] = (So * (v[i-1] * Math.cos(theta[i-1])));
			
			accx = (((-1) * (adragx)) - (amagnusx[i-1]));
		    accy = (((-1) * (adragy)) + (amagnusy[i-1]) - g);
			
			vx = ((v[i-1] * Math.cos(theta[i-1])) + (accx * dt));
		    vy = ((v[i-1] * Math.sin(theta[i-1])) + (accy * dt));
		    v[i] = Math.sqrt((Math.pow(vx,2.0)) + (Math.pow(vy,2.0)));
		    theta[i] = Math.atan(vy/vx);
		
		    x45[i] = (2*(x[i-1])) - (x[i-2]) + (accx*(dt*dt));
			y45[i] = (2*(y[i-1])) - (y[i-2]) + (accy*(dt*dt));
	    }*/
		
	    ////////////////////////////////////////////////////////////////
		//1. Opening a file to store the x vs y data
		///////////////////////////////////////////////////////////////

		// Open log file to write the angle vs time data
  	  	String filename = "GolfSimulationOutput.txt";
		PrintWriter outputFile = null;
  	  	try
  	  	{
  	  	  outputFile = new PrintWriter(new FileOutputStream("angle_vs_time.txt",false));
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

	   		plot(thetaMax, x, y); //Method below
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// 4. Find and then print period as a function of amplitude
		//////////////////////////////////////////////////////////////////////////////
		System.out.println("x	y");
		System.out.println("(m)		(m)");
		System.out.println(xMax + "  0 ");
		
		///////////////////////////////////////////////////////////////////////
		// 5. Closing file
		///////////////////////////////////////////////////////////////////////
		outputFile.close();
		
}		
    ///////////////////////////////////////////////////////////////////////
  	// 7. Live graphing of angular postion as a function of time for amplitude = printamplitude
  	///////////////////////////////////////////////////////////////////////
  	public static void plot(double theta, double[] x, double[] y)
  	{
  	  	Plot2DPanel plot1 = new Plot2DPanel();

    	// define the legend position
    	plot1.addLegend("SOUTH");

    	// add a line plot to the PlotPanel
    	plot1.addLinePlot("GolfBall", x, y);
    	plot1.setAxisLabel(0,"x (m)");
    	plot1.setAxisLabel(1,"y (m)");
    	BaseLabel title1 = new BaseLabel("x v.s. y for max angle "
    	  + theta + " degrees", Color.BLACK, 0.5, 1.1);

    	title1.setFont(new Font("Courier", Font.BOLD, 14));
    	plot1.addPlotable(title1);

    	// put the PlotPanel in a JFrame like a JPanel
    	JFrame frame1 = new JFrame("a plot panel");
    	frame1.setSize(1200, 1200);
    	frame1.setContentPane(plot1);
    	frame1.setVisible(true);
		frame1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  	}//plot
		
		
		/*System.out.printf("The angle in degrees to launch the golf ball at for maximum distance is %.3f\n",thetaMax);
		System.out.printf("The maximum distance travelled by the ball in meters is %.3f\n",xMax);*/
}
		
		