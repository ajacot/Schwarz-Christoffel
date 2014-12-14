package map

import spire.math.Complex
import spire.implicits._
import jvx.numeric.PnMatrix;

import scala.swing._
import java.awt.{ Graphics2D, Color }

import spire.math.Complex.doubleToComplex
import spire.math.Complex.intToComplex

@SerialVersionUID(102L)
class CRDTMap (val triangulation : Triangulation, val cross_ratios : Array[Double], val automorphisms : Array[(Moebius, Moebius, Moebius, Moebius)],
		       val inner_angles : IndexedSeq[Double], val Cs : Array[Complex[Double]], val As : Array[Complex[Double]]) extends Serializable {
	
	def getDiskMap(diag : Int) : DiskMap = {
	  //val pre_vertices = new Array[Complex[Double]](inner_angles.length);
	  val pre_vertices = new Array[Complex[Double]](inner_angles.length);
	
	  triangulation.get_pre_angles(diag, cross_ratios, pre_vertices);
	  
	  return new DiskMap(inner_angles, pre_vertices, Cs(diag), As(diag));
	}
	
	def getDiskMapNoAC(diag : Int) : DiskMap = {
	  val pre_vertices = new Array[Complex[Double]](inner_angles.length);
	  triangulation.get_pre_angles(diag, cross_ratios, pre_vertices);
	  return new DiskMap(inner_angles, pre_vertices, 1.0, 0);
	}
	
	
	def condition(vertices : IndexedSeq[Complex[Double]]) : (Int) => Double = {
	  return (i_diag) => {
	    val diag = triangulation.diagonals(i_diag);
	    Math.log(Triangulation.cross_ratio(vertices(diag.a), vertices(diag.b), vertices(diag.c), vertices(diag.d)).abs);
	  }
	      
	}

	def condition_lagrangien(poly : Polygon, qual : Int = 500) : Array[Array[Double]] = {
	  val k = triangulation.diagonals.length;
	  val mat = Array.ofDim[Double](k, k);
	  
	  val eps = 0.2;
	  
	  val vertices = new Array[Complex[Double]](poly.n);
	  setConstantsOnlyC(0, vertices, qual);
	    
	  val cond0 = new Array[Double](k);
	  for(i_diag <- 0 until k){
	    	val diag = triangulation.diagonals(i_diag);
	    	cond0(i_diag) = Math.log(Triangulation.cross_ratio(vertices(diag.a), vertices(diag.b), vertices(diag.c), vertices(diag.d)).abs);
	    }
	  
	  for(i_pre <- 0 until k){
	    val last_cr = cross_ratios(i_pre);
	    cross_ratios(i_pre) = -Math.exp(Math.log(-last_cr) + eps);
	    
	    setConstantsOnlyC(0, vertices, qual);
	  
	    for(i_diag <- 0 until k){
	    	val diag = triangulation.diagonals(i_diag);
	    	mat(i_diag)(i_pre) = (Math.log(Triangulation.cross_ratio(vertices(diag.a), vertices(diag.b), vertices(diag.c), vertices(diag.d)).abs) - cond0(i_diag)) / eps;
	    }
	    
	    cross_ratios(i_pre) = last_cr;
	  }
	  
	  return mat;
	}
	
	def newton(poly : Polygon, eps : Double = 0.007, qual : Int = 50){
	  
	  val k = triangulation.diagonals.length;
	  val goal = new Array[Double](k);
	  for(id <- 0 until k) goal(id) = condition(poly.points)(id);
	  
	  val vertices0 = new Array[Complex[Double]](poly.n)
	  val cond0 = new Array[Double](k);
	  for(limit <- 0 until 200){
	    
	    setConstantsOnlyC(0, vertices0, 200);
	    cond0(0) = goal(0) - condition(vertices0)(0);
	    
	    var max = Math.abs(cond0(0));
	    for(id <- 0 until k){
	      cond0(id) = goal(id) - condition(vertices0)(id);
	      max = Math.max(max, Math.abs(cond0(id)));
	    }
	    println(max);
		  
	    if(max < eps){
		    for(id <- 0 until triangulation.diagonals.length){
		      val map = getDiskMapNoAC(id);
		      val diag = triangulation.diagonals(id);
		      Cs(id) = (poly.points(diag.c) - poly.points(diag.a)) / map.pre_apply(diag.a, diag.c);
		      As(id) = poly.points(diag.a) - map.pre_apply(diag.a) * Cs(id);
		    }
		    return;
	    }
	    
	    val lagrangien : Array[Array[Double]] = condition_lagrangien(poly, qual);//(limit+4) * 25);
	    
	    val indx = new Array[Int](k);
      
	    PnMatrix.ludcmp(lagrangien, k, indx);
	    PnMatrix.lubksb(lagrangien, k, indx, cond0);
      
	    for(id <- 0 until k) cross_ratios(id) = -Math.exp(Math.log(-cross_ratios(id)) + cond0(id));
	  
	  }
	  
	  setConstants(0, null, 200);
	  
	  val map0 = getDiskMapNoAC(0);
	  val diag0 = triangulation.diagonals(0);
	  val C0 = (poly.points(diag0.c) - poly.points(diag0.a)) / (Cs(0) *map0.pre_apply(diag0.a, diag0.c));
	  val A0 = poly.points(diag0.a) - (map0.pre_apply(diag0.a) * Cs(0) + As(0)) * C0;
    
	  for(id <- 0 until triangulation.diagonals.length){
		  Cs(id) = C0 * Cs(id);
		  As(id) = As(id) * C0 + A0;
	  }
	}
	
	def setConstants(start : Int = 0, points : Array[Complex[Double]] = null, qual : Int = 500) = triangulation.setConstants(start, cross_ratios, inner_angles, Cs, As, points, qual)
	def setConstantsOnlyC(start : Int = 0, points : Array[Complex[Double]] = null, qual : Int = 500) = triangulation.setConstantsOnlyC(start, cross_ratios, inner_angles, Cs, points, qual);
	
	def automorphism(from : Int, to : Int) : Moebius = {
	  if(from == to) return Moebius.identity;
	  
	  def next(i_diag : Int) : Option[Moebius] = {
	    if(i_diag > 0){
	      if(i_diag - 1 == to) return Some(Moebius.identity);
	      val diag = triangulation.diagonals(i_diag - 1);
	      val autos = automorphisms(i_diag - 1)
	      return (next(diag.ab).map(_.compose(autos._1))) orElse (next(diag.bc).map(_.compose(autos._2)));
	    }else if(i_diag < 0){
	      if(-i_diag - 1 == to) return Some(Moebius.identity);
	      val diag = triangulation.diagonals(-i_diag - 1);
	      val autos = automorphisms(-i_diag - 1)
	      return (next(diag.cd).map(_.compose(autos._3))) orElse (next(diag.da).map(_.compose(autos._4)));
	    }else{
	      return None;
	    }
	  }
	  
	  val diag = triangulation.diagonals(from);
	  val autos = automorphisms(from);
	  return ((next(diag.ab).map(_.compose(autos._1))) 
	   orElse (next(diag.bc).map(_.compose(autos._2))) 
	   orElse (next(diag.cd).map(_.compose(autos._3)))
	   orElse (next(diag.da).map(_.compose(autos._4)))
	   ).get
	}
  def best_diag(z : Complex[Double], start : Int) : (Complex[Double], Int) = {
	
    def next(i_diag : Int, z0 : Complex[Double]) : (Complex[Double], Int, Double) = {
      if(i_diag > 0){
        val diag = triangulation.diagonals(i_diag - 1);
	    val autos = automorphisms(i_diag - 1);
	    val bestAB = if(autos._1 != null) {val w = autos._1(z0); Stream(next(diag.ab, w))} else Stream();
	    val bestBC = if(autos._2 != null) {val w = autos._2(z0); Stream(next(diag.bc, w))} else Stream();
	    return ((z0, i_diag - 1, z0.abs) +: (bestAB ++ bestBC)).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
      }else if(i_diag < 0){
        val diag = triangulation.diagonals(-i_diag - 1);
	    val autos = automorphisms(-i_diag - 1);
	    val bestCD = if(autos._3 != null) {val w = autos._3(z0); Stream(next(diag.cd, w))} else Stream();
	    val bestDA = if(autos._4 != null) {val w = autos._4(z0); Stream(next(diag.da, w))} else Stream();
	    return ((z0, -i_diag - 1, z0.abs) +: (bestCD ++ bestDA)).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
      }else{
        return null;
      }
    }
    
    val diag = triangulation.diagonals(start);
    val autos = automorphisms(start);
    val bestAB = if(autos._1 != null) {val w = autos._1(z); Stream(next(diag.ab, w))} else Stream();
    val bestBC = if(autos._2 != null) {val w = autos._2(z); Stream(next(diag.bc, w))} else Stream();
    val bestCD = if(autos._3 != null) {val w = autos._3(z); Stream(next(diag.cd, w))} else Stream();
    val bestDA = if(autos._4 != null) {val w = autos._4(z); Stream(next(diag.da, w))} else Stream();
	val (z1, i_diag, dist) = ((z, start, z.abs) +: (bestAB ++ bestBC ++ bestCD ++ bestDA)).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
    
	return (z1, i_diag);
    /*
    def next(i_diag : Int, z0 : Complex[Double], d : Double) : (Complex[Double], Int) = {
      if(i_diag > 0){
        val diag = triangulation.diagonals(i_diag - 1);
	    val autos = automorphisms(i_diag - 1);
	    val AB = if(autos._1 != null) {val w = autos._1(z0); Stream((w, diag.ab, w.abs))} else Stream()
	    val BC = if(autos._2 != null) {val w = autos._2(z0); Stream((w, diag.bc, w.abs))} else Stream()
	    if(AB.isEmpty && BC.isEmpty) return (z0, i_diag - 1);
	    val (z1, next_diag, dist) = (AB ++ BC).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
	    
	    if(dist <= d) return next(next_diag, z1, dist);
	    else return (z0, i_diag - 1);
      }else if(i_diag < 0){
        val diag = triangulation.diagonals(-i_diag - 1);
	    val autos = automorphisms(-i_diag - 1);
	    val CD = if(autos._3 != null) {val w = autos._3(z0); Stream((w, diag.cd, w.abs))} else Stream()
	    val DA = if(autos._4 != null) {val w = autos._4(z0); Stream((w, diag.da, w.abs))} else Stream()
	    if(CD.isEmpty && DA.isEmpty) return (z0, -i_diag - 1);
	    val (z1, next_diag, dist) = (CD ++ DA).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
	    
	    if(dist <= d) return next(next_diag, z1, dist);
	    else return (z0, -i_diag - 1);
      }else{
        return (z0, -1);
      }
    }
    val diag = triangulation.diagonals(start);
    val d = z.abs;
    val autos = automorphisms(start);
    val AB = if(autos._1 != null) {val w = autos._1(z); Stream((w, diag.ab, w.abs))} else Stream()
    val BC = if(autos._2 != null) {val w = autos._2(z); Stream((w, diag.bc, w.abs))} else Stream()
    val CD = if(autos._3 != null) {val w = autos._3(z); Stream((w, diag.cd, w.abs))} else Stream()
    val DA = if(autos._4 != null) {val w = autos._4(z); Stream((w, diag.da, w.abs))} else Stream()
	if(AB.isEmpty && BC.isEmpty && CD.isEmpty && DA.isEmpty) return (z, start);
    val (z1, i_diag, dist) = (AB ++ BC ++ CD ++ DA).min(Ordering.by((t : (Complex[Double], Int, Double)) => t._3));
    if(dist <= d) return next(i_diag, z1, dist);
    else return (z, start);*/
  }
  import java.awt.Color;
  import java.awt.image.BufferedImage;
	
  def getCrowdingImage(width : Int, height : Int, from : Int, poly : Polygon) : BufferedImage = {
	val image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	
	
	val b = poly.bounds;
	val scaling = new Rect(b.x, b.y + b.height, b.width / width.toDouble, - b.height / height.toDouble);//Rect.getScaleInverse(poly.bounds, new Rect(0, 0, width.toDouble, height.toDouble));
	
	val map = getDiskMap(from);
	    
	def colorFrom(crowding : Double) : Int = {
	  val z = crowding - Math.log(map.C.abs);
	  if(z < 0){
	    new Color(255, 255 -(Math.atan(-z) / Math.PI * 256).toInt, 255 -(Math.atan(-z) / Math.PI * 256).toInt).getRGB();
	  }else{
	    new Color(255 - (Math.atan(z) / Math.PI * 256).toInt, 255 - (Math.atan(z) / Math.PI * 256).toInt, 255).getRGB();
	  }
	}
	
	for(idX <- 0 until width){
	  var last_z : Option[Complex[Double]] = None;
	  
	  for(idY <- (0 until height)){
	    val w = scaling(Complex(idX.toDouble + 0.5, idY.toDouble + 0.5));
	    val next_w = scaling(Complex(idX.toDouble + 0.5, idY.toDouble + 1.5));
	    
		val id_diag = triangulation.nearest_diag(w, poly);
		val diag = triangulation.diagonals(id_diag);
	    
	    val a = (poly.points(diag.a) - w).arg;
	    val b = (poly.points(diag.b) - w).arg;
	    val c = (poly.points(diag.c) - w).arg;
	    val d = (poly.points(diag.d) - w).arg;
	    
	    def dif_angle(x : Double, y : Double) : Double = AngleRange.toRange(y - x);
	    
	    if(dif_angle(a, b) + dif_angle(b, c) + dif_angle(c, d) > Math.PI){
	      try{
	        //val guess = map.inverse_euler_tr(w, poly, triangulation, (map.pre_vertices(diag.a) + map.pre_vertices(diag.b)) / 2.0).get;
	        //val z = map.safe_inverse_newton(w, guess).getOrElse(guess);
	        
		    //val z = map.inverse_newton(w, map.inverse_euler(w, poly, triangulation, (map.pre_vertices(diag.a) + map.pre_vertices(diag.b)) / 2.0));
		      //val z = map.inverse(w, (map.pre_vertices(diag.a) + map.pre_vertices(diag.b)) / 2.0 );
	          val z = last_z match{
		        case None => map.inverse_newton(w, map.inverse_euler_tr(w, poly, triangulation, 0.0).get);
		        case Some(guess) => map.inverse_newton(w, map.inverse_euler(w, guess));
		      }
	        
		      if(z.abs > 1.0) println(z.abs);
		      
		      last_z = if(poly.safe(w, next_w)) Some(z) else None;
		      
		      val crowding = Math.log(map.derivative(z).abs);
		      image.setRGB(idX, idY, colorFrom(crowding));
	      }catch{
	        case x : Exception => x.printStackTrace();image.setRGB(idX, idY, new Color(0, 255, 0).getRGB());
	      }
	    }else{
	      last_z = None;
	      image.setRGB(idX, idY, new Color(255, 255, 255).getRGB());
	    }
	    
	  }
	}
	
	//val diag = triangulation.nearest_diag(z, poly);
	return image;
  }
  
  def getDiskCrowdingImage(width : Int, height : Int, from : Int, poly : Polygon) : BufferedImage = {
	val image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
	
	
	val b = poly.bounds;
	val scaling = new Rect(-1, 1, 2 / width.toDouble, - 2 / height.toDouble);//Rect.getScaleInverse(poly.bounds, new Rect(0, 0, width.toDouble, height.toDouble));
	
	val map = getDiskMap(from);
	    
	def colorFrom(crowding : Double) : Int = {
	  val z = crowding - Math.log(map.C.abs);
	  if(z < 0){
	    new Color(255, 255 -(Math.atan(-z) / Math.PI * 256).toInt, 255 -(Math.atan(-z) / Math.PI * 256).toInt).getRGB();
	  }else{
	    new Color(255 - (Math.atan(z) / Math.PI * 256).toInt, 255 - (Math.atan(z) / Math.PI * 256).toInt, 255).getRGB();
	  }
	}
	
	for(idX <- 0 until width){
	  for(idY <- (0 until height)){
	    val w = scaling(Complex(idX.toDouble + 0.5, idY.toDouble + 0.5));
	    
	    if(w.abs <= 1.0){
	      val crowding = Math.log(map.derivative(w).abs);
	      image.setRGB(idX, idY, colorFrom(crowding));
	    }else{
	      image.setRGB(idX, idY, Color.WHITE.getRGB());
	    }
	    
	  }
	}
	
	return image;
  }
}

object CRDTMap{
  def fromPoly(poly : Polygon) : CRDTMap = {
    val triangulation = Triangulation.delaunay(poly);
    val cross_ratios = new Array[Double](triangulation.diagonals.length);
    for(id <- 0 until triangulation.diagonals.length){
      val a = poly.points(triangulation.diagonals(id).a);
      val b = poly.points(triangulation.diagonals(id).b);
      val c = poly.points(triangulation.diagonals(id).c);
      val d = poly.points(triangulation.diagonals(id).d);
      cross_ratios(id) = -Triangulation.cross_ratio(a, b, c, d).abs; 
    }
    
    val Cs = new Array[Complex[Double]](triangulation.diagonals.length);
    val As = new Array[Complex[Double]](triangulation.diagonals.length);
    
    val automorphisms = new Array[(Moebius, Moebius, Moebius, Moebius)](poly.n - 3);
    
    val crdt = new CRDTMap(triangulation, cross_ratios, automorphisms, poly.inner_angles.toArray, Cs, As);
    
    crdt.newton(poly);
    
    def next(i_diag : Int, a : Complex[Double], b : Complex[Double], c : Complex[Double]) : Moebius = {
      if(i_diag > 0){
        val diag = triangulation.diagonals(i_diag - 1);
        
        val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(i_diag - 1));
        val a1 = Complex.polar(1.0, idA);
        val b1 = Complex.polar(1.0, idB);
        val c1 = Complex.polar(1.0, idC);
    	val d1 = Complex.polar(1.0, idD);
    	
    	return Moebius.fromToPoints(a, b, c, c1, d1, a1);
      }else if(i_diag < 0){
        val diag = triangulation.diagonals(-i_diag - 1);
        
        val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(-i_diag - 1));
        val a1 = Complex.polar(1.0, idA);
        val b1 = Complex.polar(1.0, idB);
        val c1 = Complex.polar(1.0, idC);
    	val d1 = Complex.polar(1.0, idD);
    	
    	return Moebius.fromToPoints(a, b, c, a1, b1, c1);
      }else{
        return null;
      }
    }
    
    for(i_diag <- 0 until poly.n -3){
      val diag = triangulation.diagonals(i_diag);
      
      val (idA, idB, idC, idD) = Triangulation.from_cross_ratio(cross_ratios(i_diag));
      val a = Complex.polar(1.0, idA);
      val b = Complex.polar(1.0, idB);
      val c = Complex.polar(1.0, idC);
      val d = Complex.polar(1.0, idD);
      
      automorphisms(i_diag) = (next(diag.ab, b, c, a), next(diag.bc, c, a, b), next(diag.cd, d, a, c), next(diag.da, a, c, d));
    }
    
    return crdt;
  }
}

import scala.util.Random

object SwingCRDT extends SimpleSwingApplication {
  /*val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps);*/
  
  val poly = Polygon.labyrinth(4, 4, 1.0, new Random(128))._1;
  val crdtmap = CRDTMap.fromPoly(poly);
   
  //val map = crdtmap.getDiskMap(7);
  
  val d0 = crdtmap.getDiskMap(4).getDisplay;
  val d1 = crdtmap.getDiskMap(11).getDisplay;
  
  //val d1 = crdtmap.getDiskMap(0).getDisplay;
  
  val display = new DisplayArray(Array(d0, d1));
  //val display = new DisplayArray(Array(Display.round_grid((x) => x, 1.0, 0.0, 20, 10), 
  //    new Points(map.pre_vertices)));
  val b = poly.bounds;
  //val b = new Rect(-1, -1, 2, 2);
    
  val bord : Double = 50;
  
  //val image = crdtmap.getDiskCrowdingImage(300, 300, 7, poly);
  
  def top = new MainFrame {
    title = "hello";
    
    preferredSize = new Dimension(700, 700);
    
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
    	  //g.setColor(new Color(0, 255, 0, 40));
    	  
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  /*
    	  val img_x = scaling(Complex(b.x, b.y + b.height));
    	  val img_d = scaling(Complex(b.x + b.width, b.y)) - img_x;
    	  g.drawImage(image, img_x.real.toInt, img_x.imag.toInt, img_d.real.toInt, img_d.imag.toInt, null);
    	  */
    	  //g.setColor(new Color(0, 0, 0, 100));
    	  
    	  display.display(g, scaling);
    	  
    	  //g.setColor(new Color(0, 255, 0, 255));
    	  poly.getDisplay.display(g, scaling);
      }
    }
  }
}