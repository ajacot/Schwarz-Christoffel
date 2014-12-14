package plots

import map._;

import spire.math.Complex
import spire.implicits._

import scala.swing._
import java.awt.{ Graphics2D, Color }


object CrossComplex extends SimpleSwingApplication {
  val bord : Double = 50;
  
  //val rho = 1.0;
  
  val a = Complex[Double](-1.0, 1.0);
  val b = Complex[Double](0.0, 0.0);
  val c = Complex[Double](1.0, 0.0);
  
  val f = (rho : Complex[Double]) => {val h = rho * (b - a) / (b - c); (h * c + a) / (h + 1.0)};
  
  //val steps = 100;
  //val ps = (0 until steps) map ((id) => f(Math.PI * 2.0 * id.toDouble / steps.toDouble))
  
  val display = new DisplayArray(Array(
		  new Path(Array(a, b, c)), 
		  Display.square_grid((x) => f(x.exp), new Rect(-Math.PI, 0, 2 * Math.PI, 2 * Math.PI))//new ClosedPath(ps)
  ));
    
  val bound : Rect = new Rect(-1.5, -1.5, 3, 3);
  
  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
        
    	  val scaling = Rect.getScaleInverse(bound, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  display.display(g, scaling);
      }
    }
  }
  /*val poly = new Polygon(Array(NonInfP(Complex(0, 0)), NonInfP(Complex(1, 0)), NonInfP(Complex(1, 1)), NonInfP(Complex(2, 1)), InfP(Math.PI * 1.5, Math.PI * 0.5)));
  val bord : Double = 50;
  val display = poly.getDisplay;
  val b : Rect = display.bounds;
  
  println(List(1, 2, 3).foldRight(List.empty[Int])((x, y) => x +: y));
  
  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
        
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  display.display(g, scaling);
      }
    }
  }*/
}

