package map

import spire.math.Complex
import spire.implicits._
import scala.swing._
import scala.Option.option2Iterable
import spire.math.Complex.intToComplex

@SerialVersionUID(103L)
class Triangulation(val diagonals : IndexedSeq[Diagonal]) extends Serializable{
  
  //var displays : Seq[Display] = Seq.empty;
  
  def getDisplay(poly : Polygon) : Display = {
    val displays = new Array[Display](diagonals.length+1);
    
    for((diag, id) <- diagonals.view.zipWithIndex){
      displays(id) = new LineText(poly.points(diag.a), poly.points(diag.c), id.toString);
    }
    
    displays(diagonals.length) = poly.getDisplay();
    
    return new DisplayArray(displays);
  }
  
  // get the unique embedding (up to rotation) with the quadrilateral around the diagonal start a rectangle, and having the good cross-ratios.
  def get_pre_angles(start : Int, cross_ratios : IndexedSeq[Double], pre_vertices : Array[Complex[Double]]){
    def next(ie : Int){
      if(ie > 0){
        val diag = diagonals(ie-1);
        pre_vertices(diag.b) = Triangulation.from_cross_ratio(cross_ratios(ie-1), pre_vertices(diag.c), pre_vertices(diag.d), pre_vertices(diag.a));
        next(diag.ab);
        next(diag.bc);
      }else if(ie < 0){
        val diag = diagonals(-ie-1);
        pre_vertices(diag.d) = Triangulation.from_cross_ratio(cross_ratios(-ie-1), pre_vertices(diag.a), pre_vertices(diag.b), pre_vertices(diag.c));
        next(diag.cd);
        next(diag.da);
      }
    }
    
    
    val diag = diagonals(start);
    val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(start));
    val a = Complex.polar(1.0, idA);
    val b = Complex.polar(1.0, idB);
    val c = Complex.polar(1.0, idC);
    val d = Complex.polar(1.0, idD);
    
    
    pre_vertices(diag.a) = a;
    pre_vertices(diag.b) = b;
    pre_vertices(diag.c) = c;
    pre_vertices(diag.d) = d;
    
    next(diag.ab);
    next(diag.bc);
    next(diag.cd);
    next(diag.da);
    
    
  }
  
  def transpose(z : Complex[Double], start : Int, cross_ratios : IndexedSeq[Double]) : Stream[(Int, Complex[Double])] = {
    def next(last_z : Complex[Double], i_diag : Int, a : Complex[Double], b : Complex[Double], c : Complex[Double]) : Stream[(Int, Complex[Double])] = {
      if(i_diag > 0){
        val diag = diagonals(i_diag - 1);
        
        val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(i_diag - 1));
        val a1 = Complex.polar(1.0, idA);
        val b1 = Complex.polar(1.0, idB);
        val c1 = Complex.polar(1.0, idC);
    	val d1 = Complex.polar(1.0, idD);
    	
    	val moebius = Moebius.fromToPoints(a, b, c, c1, d1, a1);
    	val next_z = moebius(last_z);
    	
    	return (i_diag - 1, next_z) +: (next(next_z, diag.ab, b1, c1, a1) ++ next(next_z, diag.bc, c1, a1, b1));
      }else if(i_diag < 0){
        val diag = diagonals(-i_diag - 1);
        
        val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(-i_diag - 1));
        val a1 = Complex.polar(1.0, idA);
        val b1 = Complex.polar(1.0, idB);
        val c1 = Complex.polar(1.0, idC);
    	val d1 = Complex.polar(1.0, idD);
    	
    	val moebius = Moebius.fromToPoints(a, b, c, c1, b1, a1);
    	val next_z = moebius(last_z);
    	
    	return (i_diag - 1, next_z) +: (next(next_z, diag.cd, d1, a1, c1) ++ next(next_z, diag.da, a1, c1, d1));
      }else{
        return Stream.empty;
      }
    }
    
    
    val diag = diagonals(start);
    val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(start));
    val a = Complex.polar(1.0, idA);
    val b = Complex.polar(1.0, idB);
    val c = Complex.polar(1.0, idC);
    val d = Complex.polar(1.0, idD);
    
    return (start, z) +: (next(z, diag.ab, b, c, a) ++ next(z, diag.bc, c, a, b) ++ next(z, diag.cd, d, a, c) ++ next(z, diag.da, a, c, d));
  }
  
  def getCuts(poly : Polygon, idA : Int, idB : Int, lambda : Double = 2.0) : Int = Math.ceil((poly.points(idB) - poly.points(idA)).abs / (lambda*getGeodesicDistance(poly, idA, idB)._1)).toInt
  
  // returns the nearest non neighbouring edge to the edge (idA, idB).
  // Also returns the "parts" for each triangle.
  def getGeodesicDistance(poly : Polygon, idA : Int, idB : Int) : (Double, Array[Seq[PathDiv]]) = { // , Array[Double], Array[Seq[PathDiv]], Array[Int]
    //val distances = new Array[Double](poly.n);
    val info = new Array[Seq[PathDiv]](poly.n);
    //val back = new Array[Int](poly.n);
    val else_dist = -10.0;
    var min_dist : Double = -1;
    
    def leftSide(theta : Double, id : Int, x : Complex[Double], dist : Double) : Seq[PathDiv] = {
      if(poly.inner_ranges(id).isInClosed(theta)) Stream(
          new Cone(x, id, AngleRange.unsafe(theta, poly.inner_ranges(id).to), dist), 
          new One(poly.points(id), id, poly.points(poly.inRange(id-1)), dist)) 
      else Stream()
      /*val r = AngleRange.unsafe(theta, poly.inner_ranges(id).to); 
      if(r.size < Math.PI) Stream(
          new Cone(x, id, r, dist), 
          new One(poly.points(id), id, poly.points(poly.inRange(id-1)), dist)) 
      else Stream()*/
    }
    def rightSide(theta : Double, id : Int, x : Complex[Double], dist : Double) : Seq[PathDiv] = {
      if(poly.inner_ranges(id).isInClosed(theta)) Stream(
          new Cone(x, id, AngleRange.unsafe(poly.inner_ranges(id).from, theta), dist), 
          new One(poly.points(id), id, poly.points(poly.inRange(id+1)), dist)) 
      else Stream()
      /*val r = AngleRange.unsafe(poly.inner_ranges(id).from, theta); 
      if(r.size < Math.PI) Stream(
          new Cone(x, id, r, dist), 
          new One(poly.points(id), id, poly.points(poly.inRange(id+1)), dist)) 
      else Stream()*/
    }
    
    def next(i_diag : Int, idX : Int, idY : Int, divs : Seq[PathDiv]){
      if(i_diag > 0){
        val diag = diagonals(i_diag-1);
        info(i_diag-1) = divs;
        val x = poly.points(diag.b);
        val ((dist, theta), div0) = divs.flatMap((div) => div.dist(x).map((_, div))).headOption.getOrElse(((else_dist, -1.0), null));
        //distances(diag.b) = dist;
        //back(diag.b) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs.map(_.cut(diag.a, diag.b, diag.c, poly));
        val ab_divs = new_divs.flatMap(_._1);
        val bc_divs = new_divs.flatMap(_._2);
        
        next(diag.ab, diag.a, diag.b, ab_divs ++ leftSide(theta, diag.b, x, dist));
        next(diag.bc, diag.b, diag.c, bc_divs ++ rightSide(theta, diag.b, x, dist));
        
      }else if(i_diag < 0){
        val diag = diagonals(-i_diag-1);
        info(-i_diag-1) = divs;
        val x = poly.points(diag.d);
        val ((dist, theta), div0) = divs.flatMap((div) => div.dist(x).map((_, div))).headOption.getOrElse(((else_dist, -1.0), null));
        //distances(diag.d) = dist;
        //back(diag.d) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs.map(_.cut(diag.c, diag.d, diag.a, poly));
        val cd_divs = new_divs.flatMap(_._1);
        val da_divs = new_divs.flatMap(_._2);
        
        next(diag.cd, diag.c, diag.d, cd_divs ++ leftSide(theta, diag.d, x, dist));
        next(diag.da, diag.d, diag.a, da_divs ++ rightSide(theta, diag.d, x, dist));
      }else{
        val x = poly.points(idX);
        val y = poly.points(idY);
        val dists = divs.flatMap(_.dist(x, y));
        
        if(idX != idB && idY != idA && !dists.isEmpty) {
          //println((idX, idY, dists.min));
          min_dist = Math.min(min_dist, dists.min);
        }
      }
    }
    
    val a = poly.points(idA);
    val b = poly.points(idB);
    
    val dir = ((b - a) * Complex.i[Double]).arg;
    val id_left = poly.inRange(idA-1);
    val id_right = poly.inRange(idB+1);
    val strip = new Strip(a, b - a, 0);
    val divs0 : Seq[PathDiv] = strip +: (leftSide(dir, idA, a, 0) ++ rightSide(dir, idB, b, 0));
    
    //displays = displays ++ divs0.map(_.getDisplay) :+ new Cone(b, idB, AngleRange.unsafe(poly.inner_ranges(idB).from, dir), 0).getDisplay;
    
    //back(idA) = -1;
    //back(idB) = -1;
    
    for(i_diag <- 0 until poly.n-3){
      val diag = diagonals(i_diag);
      if(diag.a == idA && diag.b == idB) {
        //println(("ab", i_diag));
        val x = poly.points(diag.c);
        val ((dist, theta), div0) = divs0.flatMap((div) => div.dist(x).map((_, div))).head;
        if(min_dist == -1) min_dist = dist;
        //distances(diag.c) = dist;
        //back(diag.c) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs0.map(_.cut(idB, diag.c, idA, poly));
        val left_divs = new_divs.flatMap(_._1);
        val right_divs = new_divs.flatMap(_._2);
        
        next(-i_diag-1, diag.a, diag.c, right_divs ++ rightSide(theta, diag.c, x, dist));
      }
      if(diag.b == idA && diag.c == idB) {
        //println(("bc", i_diag));
        val x = poly.points(diag.a);
        val ((dist, theta), div0) = divs0.flatMap((div) => div.dist(x).map((_, div))).head;
        if(min_dist == -1) min_dist = dist;
        //distances(diag.a) = dist;
        //back(diag.c) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs0.map(_.cut(idB, diag.a, idA, poly));
        val left_divs = new_divs.flatMap(_._1);
        val right_divs = new_divs.flatMap(_._2);
        
        next(-i_diag-1, diag.a, diag.c, left_divs ++ leftSide(theta, diag.a, x, dist));
      }
      if(diag.c == idA && diag.d == idB) {
        //println(("cd", i_diag));
        val x = poly.points(diag.a);
        val ((dist, theta), div0) = divs0.flatMap((div) => div.dist(x).map((_, div))).head;
        if(min_dist == -1) min_dist = dist;
        //distances(diag.a) = dist;
        //back(diag.c) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs0.map(_.cut(idB, diag.a, idA, poly));
        //val left_divs = new_divs.flatMap(_._1);
        val right_divs = new_divs.flatMap(_._2);
        
        next(i_diag+1, diag.a, diag.c, right_divs ++ rightSide(theta, diag.a, x, dist));
      }
      if(diag.d == idA && diag.a == idB) {
        //println(("da", i_diag));
        val x = poly.points(diag.c);
        val ((dist, theta), div0) = divs0.flatMap((div) => div.dist(x).map((_, div))).head;
        if(min_dist == -1) min_dist = dist;
        //distances(diag.c) = dist;
        //back(diag.c) = if(div0 != null) div0.from else -2;
        
        val new_divs = divs0.map(_.cut(idB, diag.c, idA, poly));
        val left_divs = new_divs.flatMap(_._1);
        //val right_divs = new_divs.flatMap(_._2);
        
        next(i_diag+1, diag.a, diag.c, left_divs ++ leftSide(theta, diag.c, x, dist));
      }
    }
    
    return (min_dist, info);
  }
  
  
  // finds the diagonal sich that its middle point is the nearest to z.
  def nearest_diag(z : Complex[Double], poly : Polygon) : Int = {
    val dist = (i_diag : Int) => {
      val diag = diagonals(i_diag);
      if(poly.points(diag.c) != poly.points(diag.a))
    	  ((z - poly.points(diag.a)) / (poly.points(diag.c) - poly.points(diag.a)) - Complex(0.5, 0)).abs;
      else
        100;
    }
    
    return (0 until poly.n-3).toStream.min(Ordering.by(dist));
  }
  
  // gives an alternative vector of inner angles such that all triangles of the triangulation will kind of look like
  // equilateral triangles (all angles are multiple of pi / 3).
  def equi_inner_angles : Array[Double] = {
    val out = new Array[Double](diagonals.length + 3);
    for(id <- 0 until out.length) out(id) = Math.PI / 3;
    for(diag <- diagonals) {
      out(diag.a) += Math.PI / 3;
      out(diag.c) += Math.PI / 3;
    }
    
    return out;
  }
  	// finds for each diagonal i_diag the constants Css(i_diag) and Ass(i_diag) so that the functions that we get with
  	// getDiskMap will be mapped to the same polygon.
   	// additionally, if the array positions is not null, then the vertices of the now unique polygon will be put in this array.
	def setConstants(start : Int = 0, 
					cross_ratios : IndexedSeq[Double],
					angles : IndexedSeq[Double],
					Css : Array[Complex[Double]], 
					Ass : Array[Complex[Double]], 
					positions : Array[Complex[Double]] = null, 
					qual : Int = 500){
	  val n = diagonals.length + 3;
	  
	  val set_positions = positions != null;
	  
	  val pre_vertices = Array.ofDim[Complex[Double]](n-3, n);
	  
	  get_pre_angles(start, cross_ratios, pre_vertices(start));
	  
	  def next(i_diag : Int, lastC : Complex[Double], lastA : Complex[Double], last_pre_vertices : IndexedSeq[Complex[Double]], x : Complex[Double], y : Complex[Double], z : Complex[Double]){
	    if(i_diag > 0){
	      val id = i_diag - 1;
	      val diag = diagonals(id);
	      get_pre_angles(id, cross_ratios, pre_vertices(id));
	      
	      val a = pre_vertices(id)(diag.a);
	      val b = pre_vertices(id)(diag.b);
	      val c = pre_vertices(id)(diag.c);
	      val d = pre_vertices(id)(diag.d);
	      
	      val back = Moebius.fromToPoints(c, d, a, x, y, z);
	      
	      val (constC, constA) = getConstants(back, pre_vertices(id), last_pre_vertices, angles, qual);
	      
	      Css(id) = constC * lastC;
	      Ass(id) = constA * lastC + lastA;
	      
	      if(set_positions){
	        val map = new DiskMap(angles, pre_vertices(id), Css(id), Ass(id));
	        
	        positions(diag.b) = map.apply_qual(diag.b, qual);
	      }
	      
	      next(diag.ab, Css(id), Ass(id), pre_vertices(id), b, c, a);
	      next(diag.bc, Css(id), Ass(id), pre_vertices(id), c, a, b);
	    
	  }else if(i_diag < 0){
		  val id = -i_diag - 1;
	      val diag = diagonals(id);
	      get_pre_angles(id, cross_ratios, pre_vertices(id));
	      
	      val a = pre_vertices(id)(diag.a);
	      val b = pre_vertices(id)(diag.b);
	      val c = pre_vertices(id)(diag.c);
	      val d = pre_vertices(id)(diag.d);
	      
	      val back = Moebius.fromToPoints(a, b, c, x, y, z);
	      
	      val (constC, constA) = getConstants(back, pre_vertices(id), last_pre_vertices, angles, qual);
	      
	      Css(id) = constC * lastC;
	      Ass(id) = constA * lastC + lastA;
	      
	      if(set_positions){
	        val map = new DiskMap(angles, pre_vertices(id), Css(id), Ass(id));
	        
	        positions(diag.d) = map.apply_qual(diag.d, qual);
	      }
	      
	      next(diag.cd, Css(id), Ass(id), pre_vertices(id), d, a, c);
	      next(diag.da, Css(id), Ass(id), pre_vertices(id), a, c, d);
	    }
	  }
	  val diag = diagonals(start);
      val a = pre_vertices(start)(diag.a);
      val b = pre_vertices(start)(diag.b);
      val c = pre_vertices(start)(diag.c);
      val d = pre_vertices(start)(diag.d);
      
      Css(start) = 1.0;
      Ass(start) = 0.0;
      
      if(set_positions){
        val map = new DiskMap(angles, pre_vertices(start), 1.0, 0.0);
        
        positions(diag.a) = map.apply_qual(diag.a, qual);
        positions(diag.b) = map.apply_qual(diag.b, qual);
        positions(diag.c) = map.apply_qual(diag.c, qual);
        positions(diag.d) = map.apply_qual(diag.d, qual);
      }
      
      next(diag.ab, 1.0, 0.0, pre_vertices(start), b, c, a);
      next(diag.bc, 1.0, 0.0, pre_vertices(start), c, a, b);
      next(diag.cd, 1.0, 0.0, pre_vertices(start), d, a, c);
      next(diag.da, 1.0, 0.0, pre_vertices(start), a, c, d);
	}
	
	def getConstants(m : Moebius, 
					pre_vertices1 : IndexedSeq[Complex[Double]], 
					pre_vertices2 : IndexedSeq[Complex[Double]], 
					angles : IndexedSeq[Double], 
					qual : Int) : (Complex[Double], Complex[Double]) = {
	  val C = (m.a * m.d - m.c * m.b) / (m.d * m.d) * 
			  (angles zip pre_vertices2)
			  	.map({case (alpha, z) => (1.0 - m.b / (z * m.d))**(alpha / Math.PI - 1.0)})
			  	.fold[Complex[Double]](1)((x, y) => x * y);
	  
	  val map = new DiskMap(angles, pre_vertices1, 1.0, 0.0)
	  
	  val A = -C*map.pre_apply_qual(m.inverse_apply(0), qual);
	  
	  return (C, A);
	}
	
	def setConstantsOnlyC(start : Int = 0, 
					cross_ratios : IndexedSeq[Double],
					angles : IndexedSeq[Double],
					Css : Array[Complex[Double]], 
					positions : Array[Complex[Double]] = null, 
					qual : Int = 500){
	  val n = diagonals.length + 3;
	  
	  val set_positions = positions != null;
	  
	  val pre_vertices = Array.ofDim[Complex[Double]](n-3, n);
	  
	  get_pre_angles(start, cross_ratios, pre_vertices(start));
	  
	  def next(i_diag : Int, lastC : Complex[Double], last_pre_vertices : IndexedSeq[Complex[Double]], x : Complex[Double], y : Complex[Double], z : Complex[Double]){
	    if(i_diag > 0){
	      val id = i_diag - 1;
	      val diag = diagonals(id);
	      get_pre_angles(id, cross_ratios, pre_vertices(id));
	      
	      val a = pre_vertices(id)(diag.a);
	      val b = pre_vertices(id)(diag.b);
	      val c = pre_vertices(id)(diag.c);
	      val d = pre_vertices(id)(diag.d);
	      
	      val back = Moebius.fromToPoints(c, d, a, x, y, z);
	      
	      val constC = getConstantsOnlyC(back, pre_vertices(id), last_pre_vertices, angles, qual);
	      
	      Css(id) = constC * lastC;
	      
	      if(set_positions){
	        val map = new DiskMap(angles, pre_vertices(id), Css(id), 0);
	        
	        positions(diag.b) = positions(diag.a) + Css(id) * map.pre_apply_qual(diag.a, diag.b, qual);
	      }
	      
	      next(diag.ab, Css(id), pre_vertices(id), b, c, a);
	      next(diag.bc, Css(id), pre_vertices(id), c, a, b);
	    
	  }else if(i_diag < 0){
		  val id = -i_diag - 1;
	      val diag = diagonals(id);
	      get_pre_angles(id, cross_ratios, pre_vertices(id));
	      
	      val a = pre_vertices(id)(diag.a);
	      val b = pre_vertices(id)(diag.b);
	      val c = pre_vertices(id)(diag.c);
	      val d = pre_vertices(id)(diag.d);
	      
	      val back = Moebius.fromToPoints(a, b, c, x, y, z);
	      
	      val constC = getConstantsOnlyC(back, pre_vertices(id), last_pre_vertices, angles, qual);
	      
	      Css(id) = constC * lastC;
	      
	      if(set_positions){
	        val map = new DiskMap(angles, pre_vertices(id), Css(id), 0);
	        
	        positions(diag.d) = positions(diag.a) +  Css(id) *  map.pre_apply_qual(diag.a, diag.d, qual);
	      }
	      
	      next(diag.cd, Css(id), pre_vertices(id), d, a, c);
	      next(diag.da, Css(id), pre_vertices(id), a, c, d);
	    }
	  }
	  val diag = diagonals(start);
      val a = pre_vertices(start)(diag.a);
      val b = pre_vertices(start)(diag.b);
      val c = pre_vertices(start)(diag.c);
      val d = pre_vertices(start)(diag.d);
      
      Css(start) = 1.0;
      
      if(set_positions){
        val map = new DiskMap(angles, pre_vertices(start), 1.0, 0.0);
        
        positions(diag.a) = map.apply_qual(diag.a, qual);
        positions(diag.b) = map.apply_qual(diag.b, qual);
        positions(diag.c) = map.apply_qual(diag.c, qual);
        positions(diag.d) = map.apply_qual(diag.d, qual);
      }
      
      next(diag.ab, 1.0, pre_vertices(start), b, c, a);
      next(diag.bc, 1.0, pre_vertices(start), c, a, b);
      next(diag.cd, 1.0, pre_vertices(start), d, a, c);
      next(diag.da, 1.0, pre_vertices(start), a, c, d);
	}
	
	def getConstantsOnlyC(m : Moebius, 
					pre_vertices1 : IndexedSeq[Complex[Double]], 
					pre_vertices2 : IndexedSeq[Complex[Double]], 
					angles : IndexedSeq[Double], 
					qual : Int) : Complex[Double] = {
	  val C = (m.a * m.d - m.c * m.b) / (m.d * m.d) * 
			  (angles zip pre_vertices2)
			  	.map({case (alpha, z) => (1.0 - m.b / (z * m.d))**(alpha / Math.PI - 1.0)})
			  	.fold[Complex[Double]](1)((x, y) => x * y);
	  return C;
	}
}

object Triangulation{
  
  def simple(n : Int) : Triangulation = { // a triangulation of a n-gon
    val diags : Array[Diagonal] = new Array(n - 3);
    
  	diags(0) = new Diagonal(0, 1, 2, 3, -1, -1, -1, 1);
  
    for(id <- 1 until n-4){
      diags(id) = new Diagonal(0, id+1, id+2, id+3, 
    		  				   id-1, -1, -1, id+1);
    }
    
  	diags(n-4) = new Diagonal(0, n-3, n-2, n-1, n-5, -1, -1, -1);
    
    return new Triangulation(diags);
  }
  
  def delaunay(poly : Polygon) : Triangulation = {
    val diags : Array[Diagonal] = new Array(poly.n - 3);
    
    var ie = 1;
    
  	def findNext(idA : Int, idB : Int, idE : Int) : (Int, Int, Int) = {
  	  if(idB - idA < 3) return (0, idA+1, 0);
  	  
  	  val a = poly.points(idA);
  	  val b = poly.points(idB);
  	  
  	  val safe : (Int => Boolean) = (id) => poly.inside(idA, id) && poly.inside(idB, id);
  	  
  	  val weight : (Int => Double) = (id) => {val c = poly.points(id); -Math.abs(((a - c) / (b - c)).arg);};
  	  
  	  val safes = (idA+1 until idB).filter(safe);
  	  
  	  val idC : Int = (idA+1 until idB).filter(safe).sortBy(weight).head
	  	  
  	  if(idC == idA+1){
  		  val this_ie = ie;
  		  ie+=1;
  		  val (ab, b, bc) = findNext(idA+1, idB, this_ie);
  		  diags(this_ie-1) = new Diagonal(idA+1, b, idB, idA, ab, bc, -idE, 0);
  		  
  		  return (0, idA+1, this_ie);
  	  }else if(idC == idB-1){
  		  val this_ie = ie;
  		  ie+=1;
  		  val (ab, b, bc) = findNext(idA, idB-1, this_ie);
  		  diags(this_ie-1) = new Diagonal(idA, b, idB-1, idB, ab, bc, 0, -idE);
  		  
  		  return (this_ie, idB-1, 0);
  		  
  	  }else{
	  	  
	  	  val this_ie = ie;
	  	  ie+=2;
	  	  
  		  val (ab0, b0, bc0) = findNext(idA, idC, this_ie);
  		  val (ab1, b1, bc1) = findNext(idC, idB, this_ie+1);
  		  diags(this_ie-1) = new Diagonal(idA, b0, idC, idB, ab0, bc0, this_ie+1, -idE);
  		  diags(this_ie) = new Diagonal(idC, b1, idB, idA, ab1, bc1, -idE, this_ie);
  		  
  		  return(this_ie, idC, this_ie+1);
  	  }
  	  
  	}
  
  	findNext(0, poly.n-1, 0);
  	
    return new Triangulation(diags);
  }
  
  def moebius_from_cross_ratio(cross_ratio : Double) : Moebius = 
    Moebius.fromToPoints(-1, Complex(0, -1), 1, -1, Complex.polar(1.0, from_cross_ratio_right(cross_ratio)), 1);
  
  // the angles for the rectangle inscribed in the disk having a specific cross-ratio
  def from_cross_ratio(cross_ratio : Double) : (Double, Double, Double, Double) = {
    val theta = Math.atan(Math.sqrt(-cross_ratio));
    return (theta, Math.PI - theta, theta - Math.PI, -theta);
  }
  
  def from_cross_ratio_right(cross_ratio : Double) : Double = 2 * Math.atan(cross_ratio);
  
  // returns the unique point d so that the cross_ratio(x, y, z, d) = cross_ratio
  def from_cross_ratio(cross_ratio : Double, x : Complex[Double], y : Complex[Double], z : Complex[Double]) : Complex[Double] = {
    if(x == y || y == z || z == x) return x;
    //println((x, y, z));
    val h = cross_ratio * (y - x) / (z - y);
    if(h == -1) return x;
    return (h*z + x) / (h + 1);
  }
  
  def cross_ratio(a : Complex[Double], b : Complex[Double], c : Complex[Double], d : Complex[Double]) : Complex[Double] = 
    (d - a) * (b - c) / ((c - d) * (a - b));
  
  // cross ratio for points on the disk given their angle
  def cross_ratio(a : Double, b : Double, c : Double, d : Double) : Double = 
    Math.sin((d - a) / 2) * Math.sin((b - c) / 2) / (Math.sin((c - d) / 2) * Math.sin((a - b) / 2));
}


@SerialVersionUID(108L)
class Diagonal(val a : Int, val b : Int, val c : Int, val d : Int, 
			   val ab : Int, val bc : Int, val cd : Int, val da : Int) extends Serializable;


// those are the "parts" of the interior when calculating the geodesic distance
abstract class PathDiv{
  def dist(x : Complex[Double]) : Option[(Double, Double)];
  def dist(a : Complex[Double], b : Complex[Double]) : Option[Double];
  def cut(d : Complex[Double]) : (Seq[PathDiv], Seq[PathDiv]);
  def cut(a : Int, b : Int, c : Int, poly : Polygon) : (Seq[PathDiv], Seq[PathDiv]);
  def getDisplay : Display;
  def from : Int;
}

class Cone(val origin : Complex[Double], val idO : Int, val range : AngleRange, val init_dist : Double) extends PathDiv {
  def dist(x : Complex[Double]) : Option[(Double, Double)] = {
    val d = x - origin;
    val theta = d.arg;
    if(range.isInClosed(theta)) Some((d.abs + init_dist, theta)) else None;
  }
  
  def dist(a : Complex[Double], b : Complex[Double]) : Option[Double] = {
    val d = b - a;
    val d_abs = d.abs;
    //val d_norm = d / (d_abs * d_abs);
    val dx = origin - a;
    val t_best = (d.real * dx.real + d.imag * dx.imag) / (d_abs * d_abs);
    
    val best = a + d * t_best;
    val dist_best = (best - origin).abs;
    val theta_best = (best - origin).arg;
    
    val t_from = Math.max(t_best + dist_best * Math.tan(range.from - theta_best) / d_abs, 0.0); 
    val t_to = Math.min(t_best + dist_best * Math.tan(range.to - theta_best) / d_abs, 1.0);
    
    if(t_to <= t_from) return None;
    if(t_to <= t_best) return Some(init_dist + (a + t_to*d - origin).abs);
    if(t_from >= t_best) return Some(init_dist + (a + t_from*d - origin).abs);
    return Some(init_dist + (best - origin).abs);
  }
  
  def cut(d : Complex[Double]) : (Seq[Cone], Seq[Cone]) = {
    val (range0, range1) = range.cutAt((d - origin).arg);
    return (range0.map(new Cone(origin, idO, _, init_dist)), range1.map(new Cone(origin, idO, _, init_dist)));
  }
  def cut(a : Int, b : Int, c : Int, poly : Polygon) : (Seq[PathDiv], Seq[PathDiv]) = {
    val alpha = (poly.points(a) - origin).arg;
    val beta = (poly.points(b) - origin).arg;
    val gamma = (poly.points(c) - origin).arg;
    
    //if(a == 2 && )
    val range_left = {val r = AngleRange(alpha, beta); if(r.size > Math.PI) Seq() else range.intersection(r)}
    val left = if(a == idO) Seq(this) else range_left.map(new Cone(origin, idO, _, init_dist));
    val range_right = {val r = AngleRange(beta, gamma); if(r.size > Math.PI) Seq() else range.intersection(r)}
    val right = if(c == idO) Seq(this) else range_right.map(new Cone(origin, idO, _, init_dist));
    
    return (left, right);
  }
  
  override def toString : String = "Cone(" + origin.toString + ", (" + range.from.toString + ", " + range.size.toString + "))";
  
  def getDisplay() : Display = new DisplayArray(Array(
      new Line(origin, origin + Complex.polar(1.0, range.from)), new Line(origin, origin + Complex.polar(1.0, range.to)),
      new Path((0 to 40).map((id) => origin + Complex.polar(0.5, range.from + id.toDouble * range.size / 40)))
  ));
  def from = idO;
}

class Strip(val a : Complex[Double], val d : Complex[Double], val init_dist : Double) extends PathDiv {
  val width = d.abs;
  
  def dist(x : Complex[Double]) : Option[(Double, Double)] = {
    if(width == 0) return None;
    val y = (x - a) / d; 
    if(0 <= y.real && 1 >= y.real && y.imag >= 0) Some((Math.abs(y.imag) * width + init_dist, (d * Complex.i[Double]).arg)) else None;
  }
  def dist(x : Complex[Double], y : Complex[Double]) : Option[Double] = {
	val dx = (x - a) / d;
    val dy = (y - a) / d;
    val ddx = dx.imag * width;
    val ddy = dy.imag * width;
    val t_x = dx.real;
    val t_y = dy.real;
    
    val t_from = Math.max(t_x, 0);
    val t_to = Math.min(t_y, 1);
    
    
    if(t_to < t_from) return None;
    //if(ddx < 0 || ddx < 0) return None
    if(t_to == t_from) return Some(init_dist + Math.min(ddx, ddy));
    if(ddx == ddy) return Some(init_dist + ddx);
    if(ddx < ddy) {
      val t = (t_from - t_x) / (t_y - t_x);
      return Some(init_dist + ddx * (1-t) + ddy * t);
    }else{
      val t = (t_to - t_x) / (t_y - t_x);
      return Some(init_dist + ddx * (1-t) + ddy * t);
    }
  }
  def cut(x : Complex[Double]) : (Seq[Strip], Seq[Strip]) = {
    if(width == 0) return (Seq(), Seq());
    val y = (x - a) / d;
    if(y.real <= 0.0){
      return (Seq(this), Seq());
    }else if(y.real < 1.0){
      val dd = d * y.real;
      return (Seq(new Strip(a+dd, d-dd, init_dist)), Seq(new Strip(a, dd, init_dist)));
    }else{
      return (Seq(), Seq(this));
    }
  }
  def cut(a : Int, b : Int, c : Int, poly : Polygon) : (Seq[PathDiv], Seq[PathDiv]) = cut(poly.points(b));
  override def toString : String = "Strip";
  
  def getDisplay() : Display = new DisplayArray(Array(
		new Line(a, a + 0.4 * d * Complex.i[Double]), new Line(a + d, a + d + 0.4 * d * Complex.i[Double])
      , new Line(a + 0.2 * d * Complex.i[Double] , a + d + 0.2 * d * Complex.i[Double])
  ));
  def from = -1;
}

// a part that only contains one point
class One(val origin : Complex[Double], val idO : Int, val d : Complex[Double], val init_dist : Double) extends PathDiv {
  lazy val dist0 = (d - origin).abs;
  lazy val dir = (d - origin).arg;
  
  def dist(x : Complex[Double]) : Option[(Double, Double)] = if(x == d) Some((init_dist + dist0, dir)) else None;
  def dist(a : Complex[Double], b : Complex[Double]) : Option[Double] = 
    if(a == d || b == d) Some(init_dist + dist0) else None;
  def cut(x : Complex[Double]) : (Seq[PathDiv], Seq[PathDiv]) = (Seq(), Seq());
  def cut(a : Int, b : Int, c : Int, poly : Polygon) : (Seq[PathDiv], Seq[PathDiv]) = 
    if(a==idO) (Seq(this), Seq())
    else if(c==idO) (Seq(), Seq(this))
    else (Seq(), Seq())
  
  def getDisplay() : Display = new Path((0 to 4).map((id) => d + Complex.polar(0.1, 2 * Math.PI * id.toDouble / 4)));
  def from = idO;
}


object SwingTriangulation extends SimpleSwingApplication {
  val bord : Double = 50;
  /*val ps = Array[Complex[Double]](
		  new Complex(0, 0), new Complex(1, 0), new Complex(2, 1), new Complex(2, 2),
		  new Complex(1, 1), new Complex(0, 1), new Complex(-1, 1), new Complex(-1, 0)
      )*/
  /*val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps);*/
  val poly = Polygon.cut(Polygon.spiral(17));
  
  
  
  val tr = Triangulation.delaunay(poly);
  
  
  val distances = tr.getGeodesicDistance(poly, 8, 9);

  /*
  println(distances._3.view.zipWithIndex.mkString("\n"));
  
  val diag = tr.diagonals(2);
  val cones = new Cone(poly.points(1), 1, poly.inner_ranges(1), 0).cut(diag.a, diag.b, diag.c, poly)._1;
  
  val display = new DisplayArray(Array(
		tr.getDisplay(poly)
	  , new PointsText(poly.points, (0 until poly.n).map((id) => id.toString + " : " + distances._1(id).toString))
      , new DisplayArray(distances._2(4).map(_.getDisplay))
	  //, new DisplayArray(cones.map(_.getDisplay))
      ));
  */
  val display = new DisplayArray(Array(
		tr.getDisplay(poly)
	  , new PointsText(poly.points, (0 until poly.n).map((id) => id.toString))
	  //, new DisplayArray(tr.displays)
      , new DisplayArray(distances._2(8).map(_.getDisplay))
	  //, new Cone(poly.points(13), 13, poly.inner_ranges(13), 0).getDisplay
	  //, new DisplayArray(cones.map(_.getDisplay))
      ));
  
  val b : Rect = Rect.boundsOf(poly.points);

  def top = new MainFrame {
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
        
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  display.display(g, scaling);
      }
    }
  }

}


