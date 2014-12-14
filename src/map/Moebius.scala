package map

import spire.math.Complex
import spire.implicits._

import scala.swing._
import spire.math.Complex.doubleToComplex
import spire.math.Complex.intToComplex


@SerialVersionUID(107L)
class Moebius(val a : Complex[Double], val b : Complex[Double], val c : Complex[Double], val d : Complex[Double]) extends Serializable {
	def apply(w : Complex[Double]) : Complex[Double] = (a * w + b) / (c * w + d);
	lazy val determinant = a * d - b * c;
	def test = determinant != 0.0;
	def inverse : Moebius = new Moebius(d, -b, -c, a);
	def inverse_apply(w : Complex[Double]) : Complex[Double] = (d * w - b) / (-c * w + a);
	def compose(m : Moebius) : Moebius = new Moebius(a*m.a + b*m.c, a*m.b + b*m.d, c*m.a + d*m.c, c*m.b + d*m.d);
	
	def getDisplay : Display = {
	  val num_peris = 6;
	  val num_radis = 16;
	  
	  val steps_peris = 80;
	  val steps_radis = 40;
	  
	  val displays = new Array[Display](num_peris + num_radis);
	  
	  var i_dis = 0;
	  
	  for(i_r <- 1 to num_peris){
	    val r = i_r.toDouble / num_peris.toDouble;
	    val ps = new Array[Complex[Double]](steps_peris);
	    for(id <- 0 until steps_peris) 
	      ps(id) = apply(Complex.polar(r, Math.PI * 2.0 * id.toDouble / steps_peris.toDouble));
	    displays(i_dis) = new ClosedPath(ps);
	    i_dis += 1;
	  }
	  
	  for(i_r <- 0 until num_radis){
	    val v = Complex.polar(1.0, Math.PI * 2.0 * i_r.toDouble / num_radis.toDouble);
	    val ps = new Array[Complex[Double]](steps_radis);
	    for(id <- 0 until steps_radis) 
	      ps(id) = apply(v * id.toDouble / steps_radis.toDouble);
	    displays(i_dis) = new Path(ps);
	    i_dis += 1;
	  }
	  
	  return new DisplayArray(displays);
	}
}

object Moebius{
  val identity = new Moebius(1.0, 0.0, 0.0, 1.0);
  def fromPoints(z0 : Complex[Double], z1 : Complex[Double], z2 : Complex[Double]) : Moebius = 
    new Moebius(z1 - z2, -z0 * (z1 - z2), z1 - z0, -z2 * (z1 - z0));
  def fromToPoints(w0 : Complex[Double], w1 : Complex[Double], w2 : Complex[Double], z0 : Complex[Double], z1 : Complex[Double], z2 : Complex[Double]) : Moebius =
    fromPoints(z0, z1, z2).inverse compose fromPoints(w0, w1, w2);
  
  def half_ratio(a : Complex[Double], b : Complex[Double], c : Complex[Double]) : Complex[Double] = (a - b) / (c - b);
  def half_ratio(a : Double, b : Double, c : Double) : Complex[Double] = 
    Complex.polar(Math.sin((b - a) / 2) / Math.sin((b - c) / 2), (a - c) / 2);
  def test_ratio(a : Double, b : Double, c : Double){
    val z0 = half_ratio(a, b, c);
    val z1 = half_ratio(Complex.polar(1.0, a), Complex.polar(1.0, b), Complex.polar(1.0, c));
    println(z0);
    println(z1);
  }
  def test_cross_ratio(a : Double, b : Double, c : Double, d : Double){
    val z0 = Triangulation.cross_ratio(a, b, c, d);
    val z1 = Triangulation.cross_ratio(Complex.polar(1.0, a), Complex.polar(1.0, b), Complex.polar(1.0, c), Complex.polar(1.0, d));
    println(z0);
    println(z1);
  }
 
  lazy val roll : Moebius = new Moebius(Complex(1, -2), -1, -1, Complex(1, 2));//fromToPoints(1.0, Complex.i, -1.0, Complex.i, -1.0, 1.0);
  
  lazy val roll_back : Moebius = new Moebius(Complex(1, 2), 1, 1, Complex(1, -2));
}

class DiskAutomorphism(val p : Complex[Double], val theta : Double, val ordre : Boolean){
  def apply(w : Complex[Double]) : Complex[Double] = {
    val v = Complex.polar(1.0, theta);
    if(ordre) v * (w - p) / (w * p.conjugate - 1);
    else (w*v - p) / (w * v * p.conjugate - 1);
  }
  def toMoebius : Moebius = {
    val v = Complex.polar(1.0, theta);
    if(ordre) return new Moebius(v, -p * v, p.conjugate, -1);
    else return new Moebius(v, -p, p.conjugate * v, -1);
  }
  def inverse : DiskAutomorphism = new DiskAutomorphism(p, -theta, !ordre);
}

object DiskAutomorphism{
  // not so good !
  lazy val roll = new DiskAutomorphism(Complex(-1.0 / 5.0, -2.0 / 5.0), 2.0 * Math.atan(2.0), false);
}


object SwingMoebius extends SimpleSwingApplication {
  val bord : Double = 50;
  
  val m = new DiskAutomorphism(Complex(0.5, 0), Math.PI, true).toMoebius;
  val mm = new DiskAutomorphism(Complex(0.5, 0), Math.PI, true).inverse.toMoebius;
  
  val mmm = m.compose(mm)
  
  val roll = Moebius.roll;//;DiskAutomorphism.roll.toMoebius;
  
  //println((roll.a, roll.b, roll.c, roll.d));
  
  println(roll(1));
  println(roll(Complex.i));
  println(roll(-1));
  
  val display = roll.getDisplay;
  val b = new Rect(-1, -1, 2, 2);
  
  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
        
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  display.display(g, scaling);
      }
    }
  }
}
