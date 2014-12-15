package map

import spire.math.Complex
import spire.implicits._
import scala.collection.immutable.Range.Double
import scala.swing._
import java.awt.Graphics2D
import scala.collection.immutable.Stream.consWrapper

import scala.util.Random
  
@SerialVersionUID(105L)
class AngleRange(val from : Double, val size : Double, private val prob : Boolean) extends Serializable{

  def to : Double = if(prob) from + size - Math.PI * 2 else from + size;
  def isInClosed(theta : Double) : Boolean = 
    if(prob) theta >= from || theta <= from + size - Math.PI * 2
  		else theta >= from && theta <= from + size
  def isInOpen(theta : Double) : Boolean = 
    if(prob) theta > from || theta < from + size - Math.PI * 2
  		else theta > from && theta < from + size
  
  		
  def cutAt(theta : Double) : (Seq[AngleRange], Seq[AngleRange]) = {
  	if(size == 0) return (Seq(), Seq());
  	val a = AngleRange.toRange(from - theta);
  	val b = AngleRange.toRange(to - theta);
  	
  	if(a == b) return (Seq(this), Seq(this));
  	if(a < b){
  		if(0 <= a) return (Seq(this), Seq(this));
  		if(0 < b) return (Seq(AngleRange.unsafe(a + theta, theta)), Seq(AngleRange.unsafe(theta, b + theta)))
  		return (Seq(this), Seq(this));
  	}
  	val back_theta = AngleRange.toRange(theta + Math.PI);
  	if(0 < b) return (Seq(AngleRange.unsafe(a + theta, theta)), Seq(AngleRange.unsafe(theta, b + theta)));
	if(0 < a) return (Seq(this), Seq(this))
	return (Seq(AngleRange.unsafe(a + theta, theta)), Seq(AngleRange.unsafe(theta, b + theta)));
  }
  def intersection(range : AngleRange) : Seq[AngleRange] = {
    if(range.size == 0) return Seq();
    if(size == 0) return Seq();
    val h_size = size * 0.5;
    val a = AngleRange.toRange(range.from - (from + h_size));
    val b = AngleRange.toRange(range.to - (from + h_size));
    
    if(a < b) {
      val x = Math.max(a, -h_size);
      val y = Math.min(b, h_size);
      if(x < y) return Seq(AngleRange.unsafe(x + from + h_size, y + from + h_size));
      return Seq();
    }
    return (if(a < h_size) Seq(AngleRange(range.from, to)) else Seq()) ++ (if(-h_size < b) Seq(AngleRange(from, range.to)) else Seq());
  }
}

object AngleRange{
  def unsafe(_a : Double, _b : Double) : AngleRange = {
    val a = toRange(_a);
    val b = toRange(_b);
    if(a == b) return new AngleRange(a, 0, false);
    if(a < b){
      new AngleRange(a, b - a, false);
    }else{
      new AngleRange(a, b - a + 2*Math.PI, true);
    }
  }
  def apply(a : Double, b : Double) : AngleRange = 
    if(a < b){
      new AngleRange(a, b - a, false);
    }else{
      new AngleRange(a, b - a + 2*Math.PI, true);
    }
  def size(from : Double, size : Double) : AngleRange = new AngleRange(from, size, from + size >= Math.PI*2);
  
  def toRange(theta : Double) : Double = { // -pi to pi
    val mm = theta % (2*Math.PI);
    if(mm <= -Math.PI) return mm + 2*Math.PI;
    if(mm > Math.PI) return mm - 2*Math.PI;
    return mm;
  }
}


@SerialVersionUID(105L)
class Polygon (val points : IndexedSeq[Complex[Double]]) extends Serializable {
	
	val n = points.length;
	
	val inner_ranges : Array[AngleRange] = Ring.getTriples(points)
			.map((p) => AngleRange(AngleRange.toRange((p._3 - p._2).arg), AngleRange.toRange((p._1 - p._2).arg))).toArray;
	val inner_angles = inner_ranges.view.map(_.size).toArray;
	
	def sum_angles : Double = inner_angles map (Math.PI - _) sum
	def test : Boolean = sum_angles - Math.PI * 2.0 == 0.0
	
	def mod(x : Double, m : Double) : Double = if(x < 0){return x%m + m}else{return x%m}
	 
	lazy val bounds : Rect = Rect.boundsOf(points);
	
	def getDisplay() : Display = new ClosedPath(points);//new DisplayArray(IndexedSeq(new ClosedPath(points), new PointsText(points, (0 until n).view.map(_.toString))));
	
  def inside(idA : Int, idB : Int) : Boolean = {
    val a = points(idA);
    val b = points(idB);
    if(Math.abs(idA-idB) == 1) return true;
    if(!inner_ranges(idA).isInClosed((b-a).arg)) {return false;}
    if(!inner_ranges(idB).isInClosed((a-b).arg)) {return false;}
    if(idA == idB) {return true}
    if(a == b) {return false}
    for((c, d) <- Ring.getPairs(points)){
      if(c != a && c != b && d != a && d != b){
        if(cross(a, b, c, d)) {return false;}
      }
    }
    return true;
  }
	
  def sideOf(a : Complex[Double], b : Complex[Double], z : Complex[Double]) : Double = 
    Math.signum((b.real - a.real)*(z.imag - b.imag) - (b.imag - a.imag)*(z.real - b.real));
  
  def difSide(s : Double, t : Double) : Boolean = s != t || s==0;
  
  def safe(a : Complex[Double], b : Complex[Double]) : Boolean = Ring.getPairs(points).forall({case (c, d) => !cross(a, b, c, d)})
  
  def cross(a : Complex[Double], b : Complex[Double], c : Complex[Double], d : Complex[Double]) : Boolean = {
    val sideC = sideOf(a, b, c);
    val sideD = sideOf(a, b, d);
    val sideA = sideOf(c, d, a);
    val sideB = sideOf(c, d, b);
    
    if(sideA == 0 && sideB == 0 && sideC == 0 && sideD == 0){ // colinear
    	val rC = ((c - a) / (b - a)).real;
    	val rD = ((d - a) / (b - a)).real;
    	return (rC >= 0 && rC <= 1) || (rD >= 0 && rD <= 1);
    }else{
    	return difSide(sideOf(a, b, c), sideOf(a, b, d)) && difSide(sideOf(c, d, a), sideOf(c, d, b));
    }
  }
  def getGrid(num : Int, steps : Int) : (IndexedSeq[Stream[IndexedSeq[Complex[Double]]]], IndexedSeq[Stream[IndexedSeq[Complex[Double]]]]) = 
  {val (x, y) = getGridE(num, steps); return (x.map(_.map(_._1)), y.map(_.map(_._1)));}
  
  // return a grid on the interior of the polygon allong with the edges where each line starts (to be able to get a better guess for the inverse function)
  def getGridE(num : Int, steps : Int) : (IndexedSeq[Stream[(IndexedSeq[Complex[Double]], Int)]], IndexedSeq[Stream[(IndexedSeq[Complex[Double]], Int)]]) = {
    
    var maxX, minX, maxY, minY : Int = 0;
    for(id <- 1 until points.length){
      val p = points(id);
      if(p.real > points(maxX).real) maxX = id;
      if(p.real < points(minX).real) minX = id;
      
      if(p.imag > points(maxY).imag) maxY = id;
      if(p.imag < points(minY).imag) minY = id;
    }
    
    val height = points(maxX).real - points(minX).real;
    val width = points(maxY).imag - points(minY).imag;
    
    val r = Math.min(width, height) / num.toDouble;
    val d = Math.min(width, height) / steps.toDouble;
    
    
    val horizontal_lines = Double.inclusive(points(minY).imag+r, points(maxY).imag, r) map ((Y) => 
      getPairs(Ring.getTriples(0 until points.length).toStream flatMap {case (i0, i1, i2) => {
        val (z0, z1, z2) = (points(i0), points(i1), points(i2));
        if(z1.imag == Y){
          val side0 = if(z0.imag == Y) z0.real < z1.real else z0.imag > Y
          val side2 = if(z2.imag == Y) z2.real > z1.real else z2.imag > Y
          if(side0 != side2) Stream((z1.real, i1))
          else Stream();
        }else if((z1.imag > Y && z0.imag < Y) || (z1.imag < Y && z0.imag > Y)){
          Stream((z0.real + (z1.real - z0.real) * (Y - z0.imag) / (z1.imag - z0.imag), i0));
        }else{
          Stream.empty;
        }
      }} sortBy(_._1)).map((range) =>
        (Double.inclusive(range._1._1, range._2._1, d).map((X) => Complex(X, Y)) :+ Complex(range._2._1, Y), range._1._2)
      ));
    
    
    val vertical_lines = Double.inclusive(points(minX).real+r, points(maxX).real, r) map ((X) => 
      getPairs(Ring.getTriples(0 until points.length).toStream flatMap {case (i0, i1, i2) => {
        val (z0, z1, z2) = (points(i0), points(i1), points(i2));
        if(z1.real == X){
          val side0 = if(z0.real == X) z0.imag < z1.imag else z0.real > X
          val side2 = if(z2.real == X) z2.imag > z1.imag else z2.real > X
          if(side0 != side2) Stream((z1.imag, i1))
          else Stream();
        }else if((z1.real > X && z0.real < X) || (z1.real < X && z0.real > X)){
          Stream((z0.imag + (z1.imag - z0.imag) * (X - z0.real) / (z1.real - z0.real), i0));
        }else{
          Stream.empty;
        }
      }} sortBy(_._1)).map((range) => 
        (Double.inclusive(range._1._1, range._2._1, d).map((Y) => Complex(X, Y)) :+ Complex(X, range._2._1), range._1._2)
      ));
    
    return (horizontal_lines, vertical_lines);
  }
  
  def getGridDisplay(num : Int, steps : Int) : Display = {
    val (h_l, v_l) = getGrid(num, steps);
    return new DisplayArray(
        (h_l ++ v_l).map((line) => 
          new DisplayArray(line.map((part) => 
            new Path(part)
          ).toArray[Path])
        ).toArray[Display]
    );
  }
  
  // utility function
  def getPairs[A](in : Stream[A]) : Stream[(A, A)] = in match {
    case Stream.Empty => return Stream();
    case x #:: Stream.Empty => println("one left");return Stream();
    case x #:: (y #:: xs) => return (x, y) #:: getPairs(xs);
  }
  
  // the mathematical modulo
  def inRange(i : Int) : Int = {
    val j = i % n;
    if(j < 0){
      return j + n;
    }else{
      return j;
    }
  }
}

// represents a permutation on 3 numbers
class Permutation3(val at1 : Int, val at2 : Int){
  def apply[A](xs : Seq[A]) : Seq[A] = insertAt(at2, insertAt(at1, Seq(xs(0)), xs(1)), xs(2));
  
  def inverse[A](xs : Seq[A]) : Seq[A] = Seq(at2, at1, 0).foldLeft((xs, Seq.empty[A]))((xys, at) => {
    (xys._1.take(at) ++ xys._1.drop(at + 1), xys._1(at) +: xys._2)
  })._2
  
  def insertAt[A](at : Int, xs : Seq[A], x : A) : Seq[A] = xs.take(at) ++ Seq(x) ++ xs.drop(at);
}

object Polygon{

  // cuts the polygon to improve the convergence of the CRDT method.
  // as lambda grows, the number of cuts diminishes
  def cut(poly : Polygon, lambda : Double = 2.0) : Polygon = {
    val tr = Triangulation.delaunay(poly);
    val cuts = Ring.getPairs(0 until poly.n).map((t) => Math.max(tr.getCuts(poly, t._1, t._2, lambda), 1)).toArray;
    println(cuts mkString " ");
    val points = new Array[Complex[Double]](cuts.view.sum);
    
    var id0 = 0;
    for((idA, idB) <- Ring.getPairs(0 until poly.n)){
      val a = poly.points(idA);
      val d = poly.points(idB) - a;
      for(id <- 0 until cuts(idB)){
        points(id0 + id) = a + id.toDouble * d / cuts(idB).toDouble;
      }
      id0 += cuts(idB);
    }
    
    return new Polygon(points);
  }
  
  // cuts the polygon and also adds some degenerate inner angles at the corresponding indices in an alternative vector of inner angles.
  def cut(poly : Polygon, ancre : IndexedSeq[Double]) : (Polygon, Array[Double]) = {
    val tr = Triangulation.delaunay(poly);
    val cuts = Ring.getPairs(0 until poly.n).map((t) => Math.max(tr.getCuts(poly, t._1, t._2), 1)).toArray;
    println(cuts mkString " ");
    val m = cuts.view.sum;
    val points = new Array[Complex[Double]](m);
    val ancre1 = new Array[Double](m);
    
    var id0 = 0;
    for((idA, idB) <- Ring.getPairs(0 until poly.n)){
      val a = poly.points(idA);
      val d = poly.points(idB) - a;
      for(id <- 0 until cuts(idB)){
        points(id0 + id) = a + id.toDouble * d / cuts(idB).toDouble;
        ancre1(id0 + id) = Math.PI;
      }
      ancre1(id0) = ancre(idA);
      id0 += cuts(idB);
    }
    
    return (new Polygon(points), ancre1);
  }
 
 // creates a random labyrinth and an alternative vector of angles "unfolding" the labyrinth.
  def labyrinth(width : Int, height : Int, dd : Double = 0.8, gen : Random = new Random()) : (Polygon, Seq[Double]) = {
    val visited : Array[Array[Boolean]] = Array.ofDim(width, height);
    
    def visit(x : Int, y : Int, dx : Int, dy : Int, from : Boolean) : (Seq[Complex[Double]], Seq[Double], Double, Double, Boolean) = {
      
      if(x < 0 || x >= width || y < 0 || y >= height) return (Seq(), Seq(), 0, 0, false);
      if(visited(x)(y)) return (Seq(), Seq(), 0, 0, false);
      visited(x)(y) = true;
      
      val permutation = new Permutation3(gen.nextInt(2), gen.nextInt(3));
      val Seq((right_points, right_angles, right_right, right_left, right_ok)
            , (front_points, front_angles, front_right, front_left, front_ok)
            , (left_points, left_angles, left_right, left_left, left_ok)) = 
        permutation.inverse(permutation(Seq((dy, -dx), (dx, dy), (-dy, dx))).map({
        case (ddx, ddy) => visit(x + ddx, y + ddy, ddx, ddy, true);
      }));
      
      val (turns0, turns1, turns2, turns3) : (Double, Double, Double, Double) = 
        if(from){
	        if(right_ok && front_ok && left_ok) (0.5, 0.5, 0.5, 0.5)
	        else if(!right_ok && !front_ok && !left_ok) (0, -0.5, -0.5, 0)
	        else if(right_ok && !front_ok && left_ok) (0.5, 0, 0, 0.5)
	        else if(!right_ok && front_ok && left_ok) (0, 0, 0.5, 0.5)
	        else if(right_ok && front_ok && !left_ok) (0.5, 0.5, 0, 0)
	        else (0, 0, 0, 0)
        }else{
        	if(right_ok && front_ok && left_ok) (0, 0.5, 0.5, 0)
	        else if(!right_ok && !front_ok && !left_ok) (-0.5, -0.5, -0.5, -0.5)
	        else if(!right_ok && front_ok && !left_ok) (-0.5, 0, 0, -0.5)
	        else if(!right_ok && !front_ok && left_ok) (-0.5, -0.5, 0, 0)
	        else if(right_ok && !front_ok && !left_ok) (0, 0, -0.5, -0.5)
	        else (0, 0, 0, 0)
        }
      
      return (right_points ++ Seq(Complex[Double](2 * x + dx * dd + dy * dd, 2 * y + dy * dd - dx * dd)) ++ 
    		  front_points ++ Seq(Complex[Double](2 * x + dx * dd - dy * dd, 2 * y + dy * dd + dx * dd)) ++ 
    		  left_points,
    		  right_angles ++ Seq(right_left + front_right + turns1) ++ 
    		  front_angles ++ Seq(front_left + left_right + turns2) ++ 
    		  left_angles,
    		  right_right + turns0, 
    		  left_left + turns3,
    		  true
          )
    }
    val (points, angles, right, left, ok) = visit(0, 0, 1, 0, false);
    return (new Polygon((Seq(Complex[Double](-1, -1)) ++ points ++ Seq(Complex[Double](-1, 1))).toArray[Complex[Double]]),
        (Seq(right) ++ angles ++ Seq(left)).map((a) => (1 + a) * Math.PI).toArray);
  }
  
  val labi_split : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 3), Complex(2, 3), Complex(2, 1), Complex(1, 1), Complex(1, 2), Complex(1, 1), Complex(2, 1), Complex(2, 3), Complex(0, 3), Complex(0, 1),
      Complex(0, 0), Complex(2, 0), Complex(3, 0), Complex(5, 0), Complex(5, 1), Complex(3, 1), Complex(3, 3), Complex(3, 1), Complex(5, 1), Complex(5, 3),
      Complex(5, 5), Complex(3.5, 5), Complex(1.5, 5), Complex(0, 5), Complex(0, 4), Complex(2, 4), Complex(4, 4), Complex(4, 2), Complex(4, 4), Complex(2, 4), Complex(0, 4)
      ));
  val ancre_labi_split = Array(
      Math.PI * 0.5, Math.PI * 1.5, Math.PI * 1.5, Math.PI, Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI, Math.PI, Math.PI,
      Math.PI, Math.PI, Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI * 1.5, Math.PI * 1.5, Math.PI, Math.PI, Math.PI,
      Math.PI, Math.PI, Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI, Math.PI, Math.PI, Math.PI, Math.PI, Math.PI * 0.5
      );
  
  val labi : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 3), Complex(2, 3), Complex(2, 1), Complex(1, 1), Complex(1, 2), Complex(1, 1), Complex(2, 1), Complex(2, 3), Complex(0, 3),
      Complex(0, 0), Complex(5, 0), Complex(5, 1), Complex(3, 1), Complex(3, 3), Complex(3, 1), Complex(5, 1),
      Complex(5, 5), Complex(0, 5), Complex(0, 4), Complex(4, 4), Complex(4, 2), Complex(4, 4), Complex(0, 4)
      ));
  
  val _labi : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 3.1), Complex(2.1, 3.1), Complex(2.1, 0.9), Complex(0.9, 0.9), Complex(1, 2), Complex(1.1, 1.1), Complex(1.9, 1.1), Complex(1.9, 2.9), Complex(0, 2.9),
      Complex(0, 0), Complex(5, 0), Complex(5, 0.9), Complex(2.9, 0.9), Complex(3, 3), Complex(3.1, 1.1), Complex(5, 1.1),
      Complex(5, 5), Complex(0, 5), Complex(0, 4.1), Complex(4.1, 4.1), Complex(4, 2), Complex(3.9, 3.9), Complex(0, 3.9)
      ));
  val ancre_labi = Array(
      Math.PI * 0.5, Math.PI * 1.5, Math.PI * 1.5, Math.PI, Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI, Math.PI,
      Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI * 1.5, Math.PI * 1.5, Math.PI, Math.PI,
      Math.PI, Math.PI * 0.5, Math.PI * 0.5, Math.PI, Math.PI, Math.PI, Math.PI * 0.5
      );
  val long : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 0), Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0), Complex(5, 0), Complex(6, 0), 
      Complex(6, 1), Complex(5, 1), Complex(4, 1), Complex(3, 1), Complex(2, 1), Complex(1, 1), Complex(0, 1)
      ));
  val longL : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 0), Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0), Complex(5, 0), Complex(6, 0), 
      Complex(6, 1), Complex(5, 1), Complex(4, 1), Complex(3, 1), Complex(2, 1), Complex(1, 1), 
      Complex(1, 2), Complex(1, 3), Complex(1, 4), Complex(1, 5), 
      Complex(0, 5), Complex(0, 4), Complex(0, 3), Complex(0, 2), Complex(0, 1)
      ));
  def sideL(side : Double) : Polygon = new Polygon(Array[Complex[Double]](
      Complex(0, 0), Complex(1 + side, 0), Complex(1 + side, 1), 
      Complex(1, 1), Complex(1, 1 + side), Complex(0, 1 + side) 
      ));
  val ancre_longL = Array(
      Math.PI, Math.PI, Math.PI, Math.PI, Math.PI, Math.PI, Math.PI * 0.5, 
      Math.PI * 0.5, Math.PI, Math.PI, Math.PI, Math.PI, Math.PI, 
      Math.PI, Math.PI, Math.PI, Math.PI * 0.5, 
      Math.PI * 0.5, Math.PI, Math.PI, Math.PI, Math.PI
      );
  val knot = new Polygon(Array[Complex[Double]](
      Complex(-1, 0), Complex(1, 0), Complex(1, 1), Complex(0, 1), Complex(0, -1), 
      Complex(0.3, -1), Complex(0.3, 0.7), Complex(0.7, 0.7), Complex(0.7, 0.3), Complex(-1, 0.3)
      ));
  def spiral(n : Int) : Polygon = new Polygon(((0 until n).map((id) => Complex(0.0, 1).pow(id + 1)*(id + 2)) 
		  								   ++ (1 to n).map((id) => Complex(0.0, 1).pow(n - id + 1)*(n - id + 1 - 1.5))).toArray[Complex[Double]]);
  def ancre_spiral(n : Int) = (0 until 2*n).map((id) => if(id == 0 || id == 1 || id == n || id == n+1) Math.PI * 0.5 else Math.PI);
}


object SwingPoly extends SimpleSwingApplication {
  val bord : Double = 50;
  /*val ps = Array[Complex[Double]](
		  new Complex(0, 0), new Complex(1, 0), new Complex(2, 1), new Complex(2, 2),
		  new Complex(1, 1), new Complex(0, 1), new Complex(-1, 1), new Complex(-1, 0)
      )*/
  /*val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps);*/
  val (poly, ancre) = Polygon.labyrinth(9, 9, 0.8, new Random(59));/*Polygon.cut(new Polygon(Array[Complex[Double]](
		  Complex(0, 0), Complex(2, 0), 
		  Complex(2, 1), Complex(0, 1)
  )));*/
  
  val tr = Triangulation.delaunay(poly);
  
  
  /*val x = Complex.polar(1.0, 0.2);
  val y = Complex.polar(1.0, 1.5);
  val z = Complex.polar(1.0, -0.5);
  val d = Triangulation.from_cross_ratio(-2, x, y, z, 3);
  println(Triangulation.cross_ratio(x, y, z, d));
  */
  
  /*
  for((diag, id) <- tr.diagonals.view.zipWithIndex){
    println(id);
    println((diag.a, diag.b, diag.c, diag.d));
    println((diag.ab, diag.bc, diag.cd, diag.da));
  }*/
  
  /*
  val cross_ratios = new Array[Double](tr.diagonals.length);
  for(id <- 0 until tr.diagonals.length){
    val diag = tr.diagonals(id);
    cross_ratios(id) = -Triangulation.cross_ratio(poly.points(diag.a), poly.points(diag.b), 
    											  poly.points(diag.c), poly.points(diag.d)).abs;
  }
  
  val pre_vertices = new Array[Complex[Double]](poly.n);
  
  tr.get_pre_angles(0, cross_ratios, pre_vertices);
  
  val labels = new Array[String](poly.n);
  for(id <- 0 until poly.n) labels(id) = id.toString;
  *//*
  val display = new PointsText(pre_vertices, labels);
  val b = Rect.boundsOf(pre_vertices);
  */
  
  val display = poly.getDisplay;//tr.getDisplay(poly);
    
  //new DisplayArray(Array(poly.getDisplay(), poly.getGridDisplay(12, 20)));
  
  val b : Rect = Rect.boundsOf(poly.points);
  
  /*
  for((diag, id) <- tr.diagonals.view.zipWithIndex){
    println(id);
    println((diag.a, diag.b, diag.c, diag.d));
    println((diag.ab, diag.bc, diag.cd, diag.da));
  }*/
  
  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
        
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


object Ring{
  def getPairs[A](seq : IndexedSeq[A]) : IndexedSeq[(A, A)] = slide(seq) zip seq;
  def getTriples[A](seq : IndexedSeq[A]) : IndexedSeq[(A, A, A)] = 
		  ((slide(seq, 1) zip seq) zip slide(seq, -1)) map((x)=>(x._1._1, x._1._2, x._2));
  def slide[A](seq : IndexedSeq[A]) : IndexedSeq[A] = seq.last +: seq.dropRight(1);
  def slide[A](seq : IndexedSeq[A], at : Int) : IndexedSeq[A] = {
	  if(at <= 0){
		  val (seq0, seq1) = seq.splitAt(-at);
		  return seq1 ++ seq0;
	  }else{
		  val (seq0, seq1) = seq.splitAt(seq.length - at);
		  return seq1 ++ seq0;
	  }
  }
}