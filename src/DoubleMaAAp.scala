package map

import spire.math.Complex
import spire.implicits._
import scala.swing._
import map.CRDTMap
import scala.swing._
import map.DiskMap
import scala.Array.fallbackCanBuildFrom
import spire.math.Complex.doubleToComplex
import spire.math.Complex.intToComplex


class DoubleMap(val to : CRDTMap, val domain : Polygon, 
				val inner_angles0 : IndexedSeq[Double], 
				val Cs0 : Array[Complex[Double]], val As0 : Array[Complex[Double]]){
  /*
  lazy val domain : Polygon = {
    val n = inner_angles0.length;
    val points = new Array[Complex[Double]](n);
    
    for(i_diag <- 0 until n-3){
      val diag = to.triangulation.diagonals(i_diag);
      
      val pre_vertices = new Array[Complex[Double]](inner_angles0.length);
      to.triangulation.get_pre_angles(i_diag, to.cross_ratios, pre_vertices);
      
      val map = new DiskMap(inner_angles0, pre_vertices, Cs0(i_diag), As0(i_diag));
      
      val tempA = points(diag.a);
      val tempC = points(diag.c);
      points(diag.a) = map.apply(diag.a);
      points(diag.c) = map.apply(diag.c);
      println((tempA, points(diag.a)));
      println((tempC, points(diag.c)));
    }
    
    new Polygon(points);
  }*/
  
  def getDiskMaps(diag : Int) : (DiskMap, DiskMap) = {
    val to_map = to.getDiskMap(diag);
    val from_map = new DiskMap(inner_angles0, to_map.pre_vertices, Cs0(diag), As0(diag));
    
    return (from_map, to_map);
  }
  
  def apply(z : Complex[Double]) : Complex[Double] = {
    val from_diag = nearest_diag(z);
    
    val pre_vertices = new Array[Complex[Double]](inner_angles0.length);
    to.triangulation.get_pre_angles(from_diag, to.cross_ratios, pre_vertices);
    val from_map = new DiskMap(inner_angles0, pre_vertices, Cs0(from_diag), As0(from_diag));
    
    val w = from_map.inverse(z);
    
    val to_map = new DiskMap(to.inner_angles, pre_vertices, to.Cs(from_diag), to.As(from_diag));
    
    return to_map(w);
  }
  
  def map_apply(zs : Array[Complex[Double]]) {
    println("map_apply");
    val pre_vertices = Array.ofDim[Complex[Double]](inner_angles0.length-3, inner_angles0.length);
    
    for(id <- 0 until inner_angles0.length-3){
      to.triangulation.get_pre_angles(id, to.cross_ratios, pre_vertices(id));
    }
    
    var i_diag = nearest_diag(zs(0));
    
    var from_map = new DiskMap(inner_angles0, pre_vertices(i_diag), Cs0(i_diag), As0(i_diag));
    var to_map = new DiskMap(to.inner_angles, pre_vertices(i_diag), to.Cs(i_diag), to.As(i_diag));
    var last_w = from_map.inverse(zs(0));
    zs(0) = to_map.apply(last_w);
    
    if(last_w.abs < 1.0){
	    val tmp = to.best_diag(last_w, i_diag);
	    if(i_diag != tmp._2){
	    	println((tmp._2, i_diag));
	    	last_w = tmp._1;
	    	i_diag = tmp._2;
	    	
	    	from_map = new DiskMap(inner_angles0, pre_vertices(i_diag), Cs0(i_diag), As0(i_diag));
	    	to_map = new DiskMap(to.inner_angles, pre_vertices(i_diag), to.Cs(i_diag), to.As(i_diag));
	    	
	    }
    }
    
    for(id <- 1 until zs.length){
      val w = from_map.inverse(zs(id), last_w);//inverse_newton(zs(id), last_w);
      if((w - last_w).abs > 0.5) println((last_w, w));
      zs(id) = zs(id-1) + to_map.apply_qual(last_w, w, 10);
      
      last_w = w;
      
      if(last_w.abs < 1.0){
	      val tmp2 = to.best_diag(last_w, i_diag);
	      if(i_diag != tmp2._2){
	    	println((tmp2._2, i_diag));
	    	last_w = tmp2._1;
	    	i_diag = tmp2._2;
	    	
	    	from_map = new DiskMap(inner_angles0, pre_vertices(i_diag), Cs0(i_diag), As0(i_diag));
	    	to_map = new DiskMap(to.inner_angles, pre_vertices(i_diag), to.Cs(i_diag), to.As(i_diag));
	    	
	      }
      }
    }
    println("end_map");
  }
  
  def getDisplay() : Display = {
    val (h_l, v_l) = domain.getGrid(3, 200);
    return new DisplayArray(
        (h_l ++ v_l).map((line) => 
          new DisplayArray(line.map((part) => {
            val part_array = part.toArray;
            map_apply(part_array);
            new Path(part_array)
          }).toArray[Path])
        ).toArray[Display]
    );
  }
  
  def nearest_diag(z : Complex[Double]) : Int = {
    val dist = (i_diag : Int) => {
      val diag = to.triangulation.diagonals(i_diag);
      ((z - domain.points(diag.a)) / (domain.points(diag.c) - domain.points(diag.a)) - Complex(0.5, 0)).abs;
    }
    
    return (0 until domain.n-3).toStream.min(Ordering.by(dist));
  }
}

object DoubleMap{
  def fromCRDT(to : CRDTMap, inner_angles0 : IndexedSeq[Double]) : DoubleMap = fromCRDT(to, inner_angles0, 0, 1, 0.0, 1.0)
  def fromCRDT(to : CRDTMap, inner_angles0 : IndexedSeq[Double], id0 : Int, id1 : Int, p0 : Complex[Double], p1 : Complex[Double]) : DoubleMap = {
    val n = inner_angles0.length;
    val points = new Array[Complex[Double]](n);
    
    val Cs0 = new Array[Complex[Double]](n-3);
    val As0 = new Array[Complex[Double]](n-3);
    
    val pre_vertices = new Array[Complex[Double]](inner_angles0.length);
      
    def next(i_diag : Int, a0 : Complex[Double], c0 : Complex[Double]){
      val qual = 1000;
      if(i_diag > 0){
    	  val diag = to.triangulation.diagonals(i_diag - 1);
    	  to.triangulation.get_pre_angles(i_diag - 1, to.cross_ratios, pre_vertices);
    	  val map = new DiskMap(inner_angles0, pre_vertices, 1.0, 0.0);
    	  
    	  val a = map.pre_apply_qual(diag.a, qual);
    	  val b = map.pre_apply_qual(diag.b, qual);
    	  val c = map.pre_apply_qual(diag.c, qual);
    	  
    	  Cs0(i_diag - 1) = (c0 - a0) / (c - a);
    	  As0(i_diag - 1) = a0 - a * Cs0(i_diag - 1);
    	  
    	  points(diag.b) = b*Cs0(i_diag - 1) + As0(i_diag - 1);
    	  
    	  next(diag.ab, points(diag.a), points(diag.b));
    	  next(diag.bc, points(diag.b), points(diag.c));
    	  
      }else if(i_diag < 0){
    	  val diag = to.triangulation.diagonals(-i_diag - 1);
    	  to.triangulation.get_pre_angles(-i_diag - 1, to.cross_ratios, pre_vertices);
    	  val map = new DiskMap(inner_angles0, pre_vertices, 1.0, 0.0);
    	  
    	  val a = map.pre_apply_qual(diag.a, qual);
    	  val c = map.pre_apply_qual(diag.c, qual);
    	  val d = map.pre_apply_qual(diag.d, qual);
    	  
    	  Cs0(-i_diag - 1) = (c0 - a0) / (a - c);
    	  As0(-i_diag - 1) = a0 - c * Cs0(-i_diag - 1);
    	  
    	  points(diag.d) = d*Cs0(-i_diag - 1) + As0(-i_diag - 1); 
    	  
    	  next(diag.cd, points(diag.c), points(diag.d));
    	  next(diag.da, points(diag.d), points(diag.a));
    	  
      }
    }
    
    val diag = to.triangulation.diagonals(0);
    to.triangulation.get_pre_angles(0, to.cross_ratios, pre_vertices);
    val map = new DiskMap(inner_angles0, pre_vertices, 1.0, 0.0);
    //val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(to.cross_ratios(0));
    val a = map.pre_apply(diag.a);
    val b = map.pre_apply(diag.b);
    val c = map.pre_apply(diag.c);
    val d = map.pre_apply(diag.d);
    
    Cs0(0) = 1.0 / (b - a);
    As0(0) = -a * Cs0(0);
    
    points(diag.a) = a*Cs0(0) + As0(0);
    points(diag.b) = b*Cs0(0) + As0(0);
    points(diag.c) = c*Cs0(0) + As0(0);
    points(diag.d) = d*Cs0(0) + As0(0);
    
    next(diag.ab, points(diag.a), points(diag.b));
    next(diag.bc, points(diag.b), points(diag.c));
    next(diag.cd, points(diag.c), points(diag.d));
    next(diag.da, points(diag.d), points(diag.a));
    
    val C = (p1 - p0) / (points(id1) - points(id0));
    val A = p0 - points(id0)*C;
    
    for(id <- 0 until n-3){
      Cs0(id) = Cs0(id)*C;
      As0(id) = As0(id)*C + A;
    }
    
    for(id <- 0 until n) points(id) = points(id) * C + A;
    
    return new DoubleMap(to, new Polygon(points), inner_angles0, Cs0, As0);
  }
}

object MainDoubleMap extends SimpleSwingApplication {
  /*val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps.map(new NonInfP(_)));*/
  val (poly, inner_angles) = Polygon.cut(Polygon.spiral(8), Polygon.ancre_spiral(8));
  val map = CRDTMap.fromPoly(poly);
  
  //val inner_angles = Polygon.ancre_labi_split;
  
  println(poly.n);
  println(inner_angles.length);
  
  val mmap = DoubleMap.fromCRDT(map, inner_angles);
  
  //val (from, to) = mmap.getDiskMaps(3);
  
  val display = mmap.getDisplay;
  val b = poly.bounds;
  //println(b);
  
  val bord : Double = 50;
  
  def top = new MainFrame {
    title = "hello";
    
    preferredSize = new Dimension(700, 700);
    
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
    	  g.setColor(new Color(0, 0, 0, 100));
    	  
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  new DisplayArray(Array(display, poly.getDisplay)).display(g, scaling);
      }
    }
  }
}
	
	def inverse(z : Complex[Double]
 : Complex[Double] = from.apply(to.inverse(z));
	def map_inverse(zs : Array[Complex[Double]]){ to.map_inverse(zs);from.map_apply(zs);}
	
	def getDisplay : Display = {
	  val b = Rect.boundsOf((0 until from.pre_angles.length).map(from.apply(_)));
	  
	  val steps : Int = 15;
	  val (incr : Double, num : Int) = if(b.width > b.height){
	    val incr = b.width / steps.toDouble;
	    (incr, steps + Math.floor(b.height / incr).toInt + 2)
val incr = b.height / steps.toDouble;
	    (incr, steps + Math.floor(b.width / incr).toInt + 2)
	  }
	  
	  val displays = new Array[Display](num);
	  val steps2 = 170;
	  
	  
	  var id = 0;
	  for(X <- Double.inclusive(b.x+incr, b.x + b.width, incr)){
		  println(X);
	    val points = Double.inclusive(b.y, b.y + b.height, b.height / steps2.toDouble)
	    				.map(new Complex(X, _)).toArray;
	    map_apply(points);
	    displays(id) = new Path(points);
	    id += 1;
	  }
	  
	  for(Y <- Double.inclusive(b.y+incr, b.y + b.height, incr)){
		println(Y);
	    val points = Double.inclusive(b.x, b.x + b.width, b.width / steps2.toDouble)
	    				.map(new Complex(_, Y)).toArray;
	    map_apply(points);
	    displays(id) = new Path(points);
	    id += 1;
	  }
	  
	  return new DisplayArray(displays);
	}
}

object DoubleMap{
	def fromPoly(poly : Polygon, ancre : Seq[Int]) : DoubleMap = {
	  val to = DiskMap.fromPoly(poly);
	  val angles = new Array[Double](ancre.length);
	  val pre_angles = new Array[Double](ancre.length);
	  
	  
	  for(id <- 0 until ancre.length){
	    angles(id) = poly.inner_angles(ancre(id));
	    pre_angles(id) = to.pre_angles(ancre(id));
	  }
	  
	  val from = new DiskMap(angles, pre_angles, 1, 0);
	  from.C = 1.0 / from.pre_apply(0, 1);
	  from.A = -from.pre_apply(0)*from.C;
	  
	  return new DoubleMap(from, to);
	}
	
	def inverse(map : DoubleMap) : DoubleMap = new DoubleMap(map.to, map.from);
}

object MainDoubleMap extends SimpleSwingApplication {
  val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps.map(new NonInfP(_)));
  val ancre = List(0, 1, 3, 4);
  val map = DoubleMap.fromPoly(poly, ancre);
  val display = map.getDisplay;
  
  val b = poly.bounds;
  //println(b);
  
  val bord : Double = 50;
  
  def top = new MainFrame {
    title = "hello";
    
    preferredSize = new Dimension(700, 700);
    
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
    	  g.setColor(new Color(0, 0, 0, 100));
    	  
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  new DisplayArray(Array(display, poly.getDisplay)).display(g, scaling);
      }
    }
  }
}
*/