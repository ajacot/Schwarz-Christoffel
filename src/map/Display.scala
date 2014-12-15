package map;

import spire.math.Complex
import spire.implicits._

import scala.swing._
import java.awt.{ Graphics2D, Color }

import scala.collection.immutable.Range.Double

import map.Ring


abstract class Display {
	def display(g : Graphics2D, scaling : Rect);
	def bounds : Rect;
}

class Line(val a : Complex[Double], val b : Complex[Double]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    g.drawLine((scaling.x + scaling.width * a.real).toInt, (scaling.y + scaling.height * a.imag).toInt, 
	    	   (scaling.x + scaling.width * b.real).toInt, (scaling.y + scaling.height * b.imag).toInt);

  }
  def bounds : Rect = new Rect(Math.min(a.real, b.real), Math.min(a.imag, b.imag), Math.abs(a.real - b.real), Math.abs(a.imag - b.imag))
}

class LineText(val a : Complex[Double], val b : Complex[Double], val txt : String) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    g.drawLine((scaling.x + scaling.width * a.real).toInt, (scaling.y + scaling.height * a.imag).toInt, 
	    	   (scaling.x + scaling.width * b.real).toInt, (scaling.y + scaling.height * b.imag).toInt);
    g.drawString(txt, (scaling.x + scaling.width * (a.real + b.real) * 0.5).toInt, (scaling.y + scaling.height * (a.imag + b.imag) * 0.5).toInt);

  }
  def bounds : Rect = new Rect(Math.min(a.real, b.real), Math.min(a.imag, b.imag), Math.abs(a.real - b.real), Math.abs(a.imag - b.imag))
}

class Points(val points : IndexedSeq[Complex[Double]]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    for(p <- points){
      val x = (scaling.x + scaling.width * p.real).toInt;
      val y = (scaling.y + scaling.height * p.imag).toInt;
      g.fillOval(x-4, y-4, 8, 8);
    }
  }
  def bounds() : Rect = Rect.boundsOf(points);
}

class PointsText(val points : Seq[Complex[Double]], val texts : Seq[String]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    for((p, txt) <- points zip texts){
      val x = (scaling.x + scaling.width * p.real).toInt;
      val y = (scaling.y + scaling.height * p.imag).toInt;
      g.drawOval(x-3, y-3, 6, 6);
      g.drawString(txt, x + 4, y);
    }
  }
  def bounds() : Rect = Rect.boundsOf(points);
}

class Path(val points : IndexedSeq[Complex[Double]]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    for(id <- 1 until points.length){
      g.drawLine((scaling.x + scaling.width * points(id-1).real).toInt, (scaling.y + scaling.height * points(id-1).imag).toInt, 
	    		   (scaling.x + scaling.width * points(id).real).toInt, (scaling.y + scaling.height * points(id).imag).toInt);
    }
  }
  def bounds() : Rect = Rect.boundsOf(points);
}

class ClosedPath(val points : IndexedSeq[Complex[Double]]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    for((a, b) <- Ring.getPairs(points)){
      g.drawLine((scaling.x + scaling.width * a.real).toInt, (scaling.y + scaling.height * a.imag).toInt, 
	    		   (scaling.x + scaling.width * b.real).toInt, (scaling.y + scaling.height * b.imag).toInt);
    }
  }
  def bounds() : Rect = Rect.boundsOf(points);
}

class DisplayArray(val displays : Seq[Display]) extends Display{
  def display(g : Graphics2D, scaling : Rect){
    for(d <- displays.filter(_ != null)) d.display(g, scaling);
  }
  def bounds() : Rect = displays.view(1, displays.length).filter(_ != null).foldRight(displays(0).bounds)((d, rect) => Rect.max(d.bounds, rect));
}

object Display{
  def square_grid(f : (Complex[Double] => Complex[Double]), rect : Rect = new Rect(0, 0, 1, 1), num : Int = 50, steps : Int = 500) : Display = {
    
	  val grid_width = rect.width / num.toDouble;
	  val steps_size = rect.height / steps.toDouble;  
	  
      val displays = new DisplayArray(
    		  Double.inclusive(rect.x, rect.x + rect.width, grid_width).map((X) => 
    		    new Path(Double.inclusive(rect.y, rect.y + rect.height, steps_size).map((Y) => f(Complex(X, Y)))))
    	   ++ Double.inclusive(rect.y, rect.y + rect.height, grid_width).map((Y) => 
    		    new Path(Double.inclusive(rect.x, rect.x + rect.width, steps_size).map((X) => f(Complex(X, Y)))))
      );
	  
	  return displays;
  }
  def round_grid(f : (Complex[Double] => Complex[Double]), radius : Double = 1.0, z0 : Complex[Double] = 0.0, num_radii : Int = 100, num_peris : Int = 50) : Display = {
    
      val displays = new Array[Display](num_radii + num_peris);
	  
	  val steps_radius = 300;
	  for(id <- 0 until num_radii){
	    val dir = Complex.polar(1.0, id.toDouble * 2 * Math.PI / num_radii.toDouble);
	    val steps = (0 to steps_radius) map (z0 + Complex(radius) * dir * _.toDouble / steps_radius.toDouble) toArray;
	    for(id <- 0 until steps.length) steps(id) = f(steps(id));
	    displays(id) = new Path(steps);
	  }
	  
	  val step_steps_peri = 100;
	  for(id <- 0 until num_peris){
	    val steps_peri = (id + 1) * step_steps_peri;
	    val r = radius * (id.toDouble + 1.0 + num_peris.toDouble * 0) / (num_peris.toDouble * 1);
	    val steps = (0 to steps_peri) map ((i) => z0 + Complex.polar(r, i.toDouble * 2 * Math.PI / steps_peri.toDouble)) toArray;
	    for(id <- 0 until steps.length) steps(id) = f(steps(id));
	    displays(num_radii + id) = new ClosedPath(steps);
	  }
	  
	  return new DisplayArray(displays);
  }
  
}

class Rect(val x : Double, val y : Double, val width : Double, val height : Double){
  
  def apply(z : Complex[Double]) : Complex[Double] = Complex(x + z.real * width, y + z.imag * height);
  def inverse(z : Complex[Double]) : Complex[Double] = Complex((z.real - x) / width, (z.imag - y) / height);
  
  def union(rect : Rect) : Rect = {
    val (nx, nwidth) = if(x > rect.x) (rect.x, Math.max(rect.width, width + x - rect.x)) 
    				   else  (x, Math.max(width, rect.width + rect.x - x));
    val (ny, nheight) = if(y > rect.y) (rect.y, Math.max(rect.height, height + y - rect.y)) 
    				    else  (y, Math.max(height, rect.height + rect.y - y));
    return new Rect(nx, ny, nwidth, nheight);
  }
  
  override def toString() : String = "Rect(" + x.toString + ", " + y.toString + ", " + width.toString + ", " + height.toString + ")"
}
object Rect{
  def max(a : Rect, b : Rect) : Rect = new Rect(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.max(a.width, b.width), Math.max(a.height, b.height))
  
  def boundsOf(ps : Seq[Complex[Double]]) : Rect = {
    if(ps.length==0) return new Rect(0, 0, 0, 0);
    
    var minX = ps(0).real;
    var maxX = ps(0).real;
    var minY = ps(0).imag;
    var maxY = ps(0).imag;
    
    for(id <- 1 until ps.length){
      minX = Math.min(minX, ps(id).real);
      maxX = Math.max(maxX, ps(id).real);
      minY = Math.min(minY, ps(id).imag);
      maxY = Math.max(maxY, ps(id).imag);
    }
    return new Rect(minX, minY, maxX - minX, maxY - minY);
  }
  def getScale(bounds : Rect, frame : Rect) : Rect = {
    val scaleX = frame.width / bounds.width;
    val scaleY = frame.height / bounds.height;
    
    if(Math.abs(scaleX) > Math.abs(scaleY)){
      return new Rect(frame.x - bounds.x*scaleY + (frame.width-bounds.width*scaleY)/2, frame.y - bounds.y*scaleY, scaleY, scaleY);
    }else{
      return new Rect(frame.x - bounds.x*scaleX, frame.y - bounds.y*scaleX + (frame.height-bounds.height*scaleX)/2, scaleX, scaleX);
    }
  }
  def getScaleInverse(bounds : Rect, frame : Rect) : Rect = {
    val scaleX = frame.width / bounds.width;
    val scaleY = frame.height / bounds.height;
    
    if(Math.abs(scaleX) > Math.abs(scaleY)){
      return new Rect(frame.x - bounds.x*scaleY + (frame.width-bounds.width*scaleY)/2, frame.y + frame.height + bounds.y*scaleY, scaleY, -scaleY);
    }else{
      return new Rect(frame.x - bounds.x*scaleX, frame.y + frame.height + bounds.y*scaleX - (frame.height-bounds.height*scaleX)/2, scaleX, -scaleX);
    }
  }
}
