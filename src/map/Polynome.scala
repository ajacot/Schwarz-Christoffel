package map

import spire.math.Complex
import spire.implicits._
import scala.swing._
import java.awt.Graphics2D
import map.CRDTMap
import map.DoubleMap
import spire.math.Complex.doubleToComplex

// this part is not really used

class Polynome (val quotients : IndexedSeq[Complex[Double]]) {
	def apply(z : Complex[Double]) : Complex[Double] = 
	  quotients.foldLeft[(Complex[Double], Complex[Double])]((0.0, 1.0))((st, q) => (st._1 + q * st._2, st._2 * z))._1;
	
}

class NewtonPolynome (val xs : IndexedSeq[Complex[Double]], val quotients : IndexedSeq[Complex[Double]]) {
	def apply(z : Complex[Double]) : Complex[Double] = 
	  (xs zip quotients).foldLeft[(Complex[Double], Complex[Double])]((0.0, 1.0))((st, q) => (st._1 + q._2 * st._2, st._2 * (z - q._1)))._1;
	
	def getDisplay(radius : Double = 1.0, z0 : Complex[Double] = 0.0) : Display = {
	  val num_radii = 100;
	  val num_peris = 50;
	  val displays = new Array[Display](num_radii + num_peris);
	  
	  
	  val steps_radius = 300;
	  for(id <- 0 until num_radii){
	    val dir = Complex.polar(1.0, id.toDouble * 2 * Math.PI / num_radii.toDouble);
	    val steps = (0 to steps_radius) map (z0 + Complex(radius) * dir * _.toDouble / steps_radius.toDouble) toArray;
	    for(id <- 0 until steps.length) steps(id) = apply(steps(id));
	    displays(id) = new Path(steps);
	  }
	  
	  val step_steps_peri = 100;
	  for(id <- 0 until num_peris){
	    val steps_peri = (id + 1) * step_steps_peri;
	    val r = radius * (id.toDouble + 1.0 + num_peris.toDouble * 0) / (num_peris.toDouble * 1);
	    val steps = (0 to steps_peri) map ((i) => z0 + Complex.polar(r, i.toDouble * 2 * Math.PI / steps_peri.toDouble)) toArray;
	    for(id <- 0 until steps.length) steps(id) = apply(steps(id));
	    displays(num_radii + id) = new ClosedPath(steps);
	  }
	  
	  return new DisplayArray(displays);
	}
}

object Polynome {
  
  def newton_interpolation(xs : IndexedSeq[Complex[Double]], ys : IndexedSeq[Complex[Double]]) : NewtonPolynome = {
    val n = xs.length;
    val dds = divided_differences(xs, ys);
    val quotients = new Array[Complex[Double]](n);

    
    quotients(0) = ys(0);
    for(id <- 1 until n){
      quotients(id) = dds((id-1) * (2 * n - id) / 2);
    }
    
    return new NewtonPolynome(xs, quotients);
  }
  
  def divided_differences(xs : IndexedSeq[Complex[Double]], ys : IndexedSeq[Complex[Double]]) : Array[Complex[Double]] = {
    val n = xs.length;
    val dds = new Array[Complex[Double]](n * (n - 1) / 2);
    
    def at(i : Int, d : Int) : Int = i + (d - 1) * (2 * n - d) / 2;
    // n - d - 1 + (d - 1) * (2 * n - d) / 2
    for(i <- 0 until n-1){
      dds(i) = (ys(i+1) - ys(i)) / (xs(i+1) - xs(i));
    }
    
    for(d <- 2 until n){
      for(i <- 0 until n - d){
        dds(at(i, d)) = (dds(at(i+1, d-1)) - dds(at(i, d-1))) / (xs(i+d) - xs(i));
      }
    }
    
    return dds;
  }
}

object SwingPolynome extends SimpleSwingApplication {
  val poly = Polygon.longL;
  val map = CRDTMap.fromPoly(poly);
  
  val inner_angles = Polygon.ancre_longL;
  
  println(poly.n);
  println(inner_angles.length);
  
  val mmap = DoubleMap.fromCRDT(map, inner_angles);
  
  val (from, to) = mmap.getDiskMaps(0);
  
  val steps = 20;
  val ws = (0 until steps).map((id) => Complex.polar((id.toDouble * 0.2) % 0.6 + 0.3, Math.PI * 2.0 * id.toDouble / steps.toDouble));
  val xs = ws.toArray;
  val ys = ws.toArray;
  from.map_apply(xs);
  to.map_apply(ys);
  
  //println(xs mkString ", ");
  //println(ys mkString ", ");
  
  val polynome = Polynome.newton_interpolation(xs, ys);
  
  val display = new DisplayArray(Array(polynome.getDisplay(0.7), poly.getDisplay, new Points(ys)));
  val b = poly.bounds;
  
  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
    	  val bord = 50;
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
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
