package map

import spire.math.Complex
import spire.implicits._

import spire.math.Complex.doubleToComplex
import spire.math.Complex.intToComplex

object DiskMap{
  
	var integration_count = 0;
	
	def log_integration(){
	  //DiskMap.integration_count+=1;
	  //println(DiskMap.integration_count);
	}
}

class DiskMap(val angles : IndexedSeq[Double], val pre_vertices : IndexedSeq[Complex[Double]], var C : Complex[Double], var A : Complex[Double]) {
	
  
	def derivative(z : Complex[Double]) :  Complex[Double] = pre_derivative(z) * C;
	def recip_derivative(z : Complex[Double]) :  Complex[Double] = pre_recip_derivative(z) / C;
	def pre_derivative(z : Complex[Double]) :  Complex[Double] = {
	  var out : Complex[Double] = 1;
	  
	  for(id <- 0 until angles.length){
	    out *= (Complex(1.0) - z / pre_vertices(id))**(angles(id)/Math.PI - 1.0);
	  }
	  return out;
	}
	def pre_recip_derivative(z : Complex[Double]) :  Complex[Double] = {
	  var out : Complex[Double] = 1;
	  
	  for(id <- 0 until angles.length){
	    out *= (Complex(1.0) - z / pre_vertices(id))**(1.0 - angles(id)/Math.PI);
	  }
	  return out;
	}
	def pre_derivative_without(z : Complex[Double], without : Int) : Complex[Double] = {
	  var out : Complex[Double] = 1;
	  
	  for(id <- 0 until angles.length){
	    if(id != without) {
	      val a = pre_vertices(id);
	      if(a  == z){
	        println("oups");
	        println(id);
	      }
	      out *= (Complex(1.0) - z / a)**(angles(id)/Math.PI - 1.0);
	    }
	  }
	  return out;
	}
	def apply(a : Complex[Double], b : Complex[Double]) : Complex[Double] = pre_apply(a, b) * C;
	def pre_apply(a : Complex[Double], b : Complex[Double]) : Complex[Double] = pre_apply_qual(a, b, 500);
	def apply_qual(a : Complex[Double], b : Complex[Double], qual  : Int) : Complex[Double] = pre_apply_qual(a, b, qual) * C;
	def pre_apply_qual(a : Complex[Double], b : Complex[Double], qual  : Int) : Complex[Double] = {
	  DiskMap.log_integration();
	  
	  if(a == b) return 0;
	  
	  val d : Complex[Double] = b - a;
	  
	  val steps : Int = Math.max(2, (d.abs * qual).toInt);
	  val cuts = halfRule(a, b, steps);
	  
	  if(cuts.isEmpty){
		  
		  var out : Complex[Double] = 0;
		  val d : Complex[Double] = b - a;
		  
		  val d_phi : Complex[Double] = d / steps;
		  
		  var last_d_phi = pre_derivative(a);
		  
		  for(id <- 0 until steps){
		    //val phi0 = a + d * ((id.toDouble + 0.0) / steps.toDouble);
		    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
		    val phi2 = a + d * ((id.toDouble + 1.0) / steps.toDouble);
		    
		    val next_d_phi = pre_derivative(phi2);
		    out += (last_d_phi + Complex(4.0) * pre_derivative(phi1) + next_d_phi) * d_phi / 6.0;
		    last_d_phi = next_d_phi;
		  }
		  
		  return out;
	  }else{
	    var last_id = cuts(0);
	    var out : Complex[Double] = pre_apply_qual(a, last_id, qual);
	  	
	    for(id <- cuts.tail){
	      out += pre_apply_qual(last_id, id, qual);
	      
	      last_id = id;
	    }
	    
	    out += pre_apply_qual(last_id, b, qual);
	    
	    return out;
	  }
	}
	
	def apply(idA : Int, idB : Int) : Complex[Double] = pre_apply(idA, idB) * C + A;
	def pre_apply(idA : Int, idB : Int) : Complex[Double] = pre_apply_qual(idA, idB, 500);
	def pre_apply_qual(idA : Int, idB : Int, qual : Int) : Complex[Double] = {
	  DiskMap.log_integration();
	  
	  
	  if(idA == idB) return 0;
	  
	  val a = pre_vertices(idA);
	  val b = pre_vertices(idB);
	  
	  var out : Complex[Double] = 0;
	  val d : Complex[Double] = b - a;
	  
	  val steps : Int = Math.max(2, (d.abs * qual).toInt);
	  val d_phi : Complex[Double] = d / steps;
	  
	  out += {
		  val z1 = a + d_phi;
		  val f0 = pre_derivative_without(a, idA);
		  val f1 = pre_derivative_without(z1, idA);
		  val r = angles(idA)/Math.PI;
		  -a * ((Complex(1.0) - z1 / a)**r) * f1 / r - a * a * (f1 - f0) * ((Complex(1.0) - z1 / a)**(r+1)) / (d_phi * r * (r+1));
	  }
	  // f'' = (1 - z / a)**(r - 1)
	  // f' = - a * (1 - z / a)**r / r
	  // f = a * a * (1 - z / a)**(r + 1) / (r* (r+1))
	  
	  for(id <- 1 until steps-1){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  
	  out += {
		  val z0 = b - d_phi;
		  val f0 = pre_derivative_without(z0, idB);
		  val f1 = pre_derivative_without(b, idB);
		  val r = angles(idB)/Math.PI;
		  b * ((Complex(1.0) - z0 / b)**r) * f0 / r + b * b * (f1 - f0) * ((Complex(1.0) - z0 / b)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  return out;
	}
	
	def apply(id : Int) : Complex[Double] = pre_apply(id) * C + A;
	def apply_qual(id : Int, qual : Int) : Complex[Double] = pre_apply_qual(id, qual) * C + A;
	def pre_apply(id : Int) : Complex[Double] = pre_apply(0.0, id);
	def pre_apply_qual(id : Int, qual : Int) : Complex[Double] = pre_apply_qual(0.0, id, qual);
	def apply(a : Complex[Double], idB : Int) : Complex[Double] = pre_apply(a, idB) * C;
	def apply(idA : Int, b : Complex[Double]) : Complex[Double] = pre_apply(idA, b) * C;
	def pre_apply(a : Complex[Double], idB : Int) : Complex[Double] = pre_apply_qual(a, idB, 500);
	def pre_apply_qual(a : Complex[Double], idB : Int, qual : Int) : Complex[Double] = {
	  DiskMap.log_integration();
	  
	  val b = pre_vertices(idB);
	  val d = b - a;
	  
	  var out : Complex[Double] = 0;
	  
	  val steps : Int = Math.max(2, (d.abs * qual).toInt);
	  val d_phi : Complex[Double] = d / steps.toDouble;
	  
	  for(id <- 0 until steps-1){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  out += {
		  val z0 = b - d_phi;
		  val f0 = pre_derivative_without(z0, idB);
		  val f1 = pre_derivative_without(b, idB);
		  val r = angles(idB)/Math.PI;
		  b * ((Complex(1.0) - z0 / b)**r) * f0 / r + b * b * (f1 - f0) * ((Complex(1.0) - z0 / b)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  return out;
	}
	
	def pre_apply(idA : Int, b : Complex[Double]) : Complex[Double] = pre_apply_qual(idA, b, 500);
	def pre_apply_qual(idA : Int, b : Complex[Double], qual : Int) : Complex[Double] = {
	  DiskMap.log_integration();
	  
	  val a = pre_vertices(idA);
	  val d = b - a;
	  
	  var out : Complex[Double] = 0;
	  
	  val steps : Int = Math.max(2, (d.abs * qual).toInt);
	  val d_phi : Complex[Double] = d / steps;
	  
	  out += {
		  val z1 = a + d_phi;
		  val f0 = pre_derivative_without(a, idA);
		  val f1 = pre_derivative_without(z1, idA);
		  val r = angles(idA)/Math.PI;
		  -a * ((Complex(1.0) - z1 / a)**r) * f1 / r - a * a * (f1 - f0) * ((Complex(1.0) - z1 / a)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  for(id <- 1 until steps){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  return out;
	}
	def apply(z : Complex[Double]) : Complex[Double] = pre_apply(z) * C + A;
	def pre_apply(z : Complex[Double]) : Complex[Double] = pre_apply_qual(z, 500);
	def pre_apply_qual(z : Complex[Double], qual : Int) : Complex[Double] = {
	  DiskMap.log_integration();
	  
	  var out : Complex[Double] = 0;
	  
	  val steps : Int = Math.max(2, qual);
	  val d_phi : Complex[Double] = z / steps;
	  
	  for(id <- 0 until steps){
	    val phi0 = z * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = z * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = z * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  return out;
	}
	
	
	def map_apply(zs : Array[Complex[Double]]){
		// here we assume a path and integrate gradually
		// replacing with the result

		var last_w = apply(zs(0));
		var last_z = zs(0);
		zs(0) = last_w;
		
		for(id <- 1 until zs.length){
			val z = zs(id);
			
			val cuts = halfRule(last_z, z, 1);
			if(cuts.length == 0){
			    val z0 = (z + last_z * 3) / 4;
			    val z1 = (z + last_z) / 2;
			    val z2 = (z * 3 + last_z) / 4;
			    last_w += (pre_derivative(z0)*2 - pre_derivative(z1) + pre_derivative(z2)*2) * C * (z - last_z) / 3;
			    
			}else{
			  var lastID = cuts.head;
			  var last = pre_vertices(cuts.head);
			  
			  if(last != last_z) last_w += apply(last_z, lastID);
			    
			  for(nextID <- cuts.tail){
			    last_w += apply(lastID, nextID);
			    lastID = nextID;
			  }
			  
			  //lastID = cuts.last;
			  last = pre_vertices(lastID);
			  if(last != z) last_w += apply(lastID, z);
			  
			}
			zs(id) = last_w;
			last_z = z;
		}
	}
	
	def halfRule(a : Complex[Double], b : Complex[Double], steps : Int) : Seq[Int] = {
	  val border = 0.5 / steps.toDouble;
	  return pre_vertices.view.map((v : Complex[Double]) => (v-a) / (b - a))
	  			.zipWithIndex.collect({ case (Complex(x, y), id)
	  				if (-border <= x && x <= 1.0+border && -border <= y && y <= border)
	  				  => (x, id)
	  			}).sortBy(_._1).map(_._2);
	}
	
	// inversing the map
	def inverse(w : Complex[Double], guess : Complex[Double] = 0) : Complex[Double] = inverse_newton(w, inverse_euler(w, guess, 300));
	def inverse_euler(w : Complex[Double], guess : Complex[Double] = 0, steps : Int = 200) : Complex[Double] = {
	   var last_z : Complex[Double] = guess;
	   var last_w = apply(last_z);
	   
	   val delta = w - last_w;
	   val h : Complex[Double] = delta / steps;
	   
	   for(id <- 0 until steps){
	     last_z += h * recip_derivative(last_z);
	   }

	   return last_z;
	}
	def inverse_newton(w : Complex[Double], guess : Complex[Double]) : Complex[Double] = 
	  safe_inverse_newton(w, guess).getOrElse(inverse_euler(w, guess, 300));
	def safe_inverse_newton(w : Complex[Double], guess : Complex[Double]) : Option[Complex[Double]] = {
	  var z = guess;
	  
	  var first_time = true;
	  
	  val err : Double = 0.0001;
	  val sqr_err = Math.pow(err, 2.0); // **2.0
	  
	  for(id <- 0 until 100){
	    if(id > 10) println(id);
	    val d = apply(z) - w;
	    if(d.norm < sqr_err) return Some(z);
	    z -= d / derivative(z);
	    
	    if(z.norm > 20.0){
	      return None;
	    }
	  }
	  
	  println("no convergence");
	  
	  return None;
	}
	def inverse_euler_tr(w : Complex[Double], poly : Polygon, tr : Triangulation, guess : Complex[Double], steps : Int = 50) : Option[Complex[Double]] = {
	  
	  val w0 = apply(guess);
	  
	  def inside(z : Complex[Double], a : Complex[Double], b : Complex[Double], c : Complex[Double]) : Boolean = 
	    ((z-a) / (b-a)).imag >= 0.0 && ((z-b) / (c-b)).imag >= 0.0 &&((z-c) / (a-c)).imag >= 0.0;
	  
	  def next(i_diag : Int) : (Option[Seq[Int]], Option[Seq[Int]]) = {
	    if(i_diag > 0){
	      val diag = tr.diagonals(i_diag - 1);
	      
	      val (ab_from, ab_to) = next(diag.ab);
	      val (bc_from, bc_to) = next(diag.bc);
	      
	      val a = poly.points(diag.a)
	      val b = poly.points(diag.b)
	      val c = poly.points(diag.c)
	      
	      val from : Option[Seq[Int]] = if(inside(w0, a, b, c)) Some(Seq()) else None;
	      val to : Option[Seq[Int]] = if(inside(w, a, b, c)) Some(Seq()) else None;
	      
	      return ((from orElse ab_from orElse bc_from).map((i_diag - 1) +: _), 
	    		  (to orElse ab_to orElse bc_to).map((i_diag - 1) +: _));
	    }else if(i_diag < 0){
	      val diag = tr.diagonals(-i_diag - 1);
	      
	      val (cd_from, cd_to) = next(diag.cd);
	      val (da_from, da_to) = next(diag.da);
	      
	      val c = poly.points(diag.c)
	      val d = poly.points(diag.d)
	      val a = poly.points(diag.a)
	      
	      val from : Option[Seq[Int]] = if(inside(w0, c, d, a)) Some(Seq()) else None;
	      val to : Option[Seq[Int]] = if(inside(w, c, d, a)) Some(Seq()) else None;
	      
	      return ((from orElse cd_from orElse da_from).map((-i_diag - 1) +: _), 
	    		  (to orElse cd_to orElse da_to).map((-i_diag - 1) +: _));
	    }
	      
	    return (None, None);
	  }
	  
	  val diag = tr.diagonals(3);
	  
      val (left_from, left_to) = next(-1);
      val (right_from, right_to) = next(1);
	  
      def merge(from : Seq[Int], to : Seq[Int]) : Seq[Int] = {
        val drop = (from zip to).takeWhile((x) => x._1 == x._2).length;
        return from.drop(drop).reverse ++ to.drop(drop);
      }
      val from = left_from orElse right_from.map(_.tail)
      val to = left_to orElse right_to.map(_.tail)
      
      if(from.isEmpty || to.isEmpty) return None;
      
      val path = merge(from.get, to.get);
      
      var last_z = guess;
      
      for(i_diag <- path){
        val diag = tr.diagonals(i_diag);
        val mid = (poly.points(diag.a) + poly.points(diag.c)) / 2.0;
        
        last_z = inverse_euler(mid, last_z, steps);
      }
      
	  return Some(inverse_euler(w, last_z, steps));
	}
	
	def map_inverse(ws : Array[Complex[Double]]){
	  if(ws.length == 0) return;
	  
	  var last_z = inverse(ws(0));
	  ws(0) = last_z;
	  for(id <- 1 until ws.length){
	    last_z = inverse_newton(ws(id), last_z);
	    ws(id) = last_z;
	  }
	}
	
	
	def getDisplay() : DisplayArray = {
	  val num_radii = 20;
	  val num_peris = 10;
	  val displays = new Array[Display](num_radii + num_peris);
	  
	  
	  val steps_radius = 300;
	  for(id <- 0 until num_radii){
	    val dir = Complex.polar(1.0, id.toDouble * 2 * Math.PI / num_radii.toDouble);
	    val steps = (0 to steps_radius) map (dir * _.toDouble / steps_radius.toDouble) toArray;
	    map_apply(steps);
	    displays(id) = new Path(steps);
	  }
	  
	  val step_steps_peri = 100;
	  for(id <- 0 until num_peris-1){
	    val steps_peri = (id + 1) * step_steps_peri;
	    val r = (id.toDouble + 1.0 + num_peris.toDouble * 0) / (num_peris.toDouble * 1);
	    val steps = (0 to steps_peri) map ((i) => Complex.polar(r, i.toDouble * 2 * Math.PI / steps_peri.toDouble)) toArray;
	    map_apply(steps);
	    displays(num_radii + id) = new ClosedPath(steps);
	  }
	  
	  return new DisplayArray(displays);
	}
	/*
	def getDisplaySmooth(d_theta : Double = 0.01) : DisplayArray = {
	  val num_radii = 100;
	  val num_peris = 50;
	  val displays = new Array[Display](num_radii + num_peris);
	  
	  
	  val steps_radius = 300;
	  for(id <- 0 until num_radii){
	    val dir = Complex.polar(1.0, id.toDouble * 2 * Math.PI / num_radii.toDouble);
	    val steps = (0 to steps_radius) map (dir * _.toDouble / steps_radius.toDouble) toArray;
	    map_smooth(steps, d_theta);
	    displays(id) = new Path(steps);
	  }
	  
	  val step_steps_peri = 100;
	  for(id <- 0 until num_peris){
	    val steps_peri = (id + 1) * step_steps_peri;
	    val r = (id.toDouble + 1.0 + num_peris.toDouble * 0) / (num_peris.toDouble * 1);
	    val steps = (0 to steps_peri) map ((i) => Complex.polar(r, i.toDouble * 2 * Math.PI / steps_peri.toDouble)) toArray;
	    map_smooth(steps, d_theta);
	    displays(num_radii + id) = new ClosedPath(steps);
	  }
	  
	  return new DisplayArray(displays);
	}
	*/
}

import scala.swing._
import java.awt.{ Graphics2D, Color }

object SwingDisk extends SimpleSwingApplication {
  /*val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps);*/
  
  //val poly = Polygon.cut(Polygon.longL);
  //val crdtmap = CRDTMap.fromPoly(poly);
  
  val ps = Array[Complex[Double]](Complex(0, 1), Complex(0, 0), Complex(1, 0), Complex(1, -2), Complex(3, -2), Complex(3, -1), Complex(2, -1)
      , Complex(2, 0), Complex(3, 0), Complex(3, 1));
  val poly = new Polygon(ps);
  val crdtmap = CRDTMap.fromPoly(poly);
  
  crdtmap.getDiskMap(2).inverse_euler_tr(Complex(0.2, 0.5), poly, crdtmap.triangulation, Complex(2.7, -1.5));
  /*
  val inner_angles = Array(Math.PI * 0.5, Math.PI * 0.5, Math.PI * 1.5, 
		  					  Math.PI * 0.5, Math.PI * 0.5, Math.PI * 0.5, Math.PI * 1.5, Math.PI * 1.5, 
		  					  Math.PI * 0.5, Math.PI * 0.5);
  
  
  val prevertices = (0 until inner_angles.length).map((id) => Complex.polar(1.0, Math.PI * 2 * id.toDouble / inner_angles.length.toDouble)).seq.toArray;
  
  val map = new DiskMap(inner_angles, 
		  				prevertices, 
		  				1, 0);
  val poly = new Polygon((0 until inner_angles.length).map((i) => map.apply(i)).seq);
  */
  
  
  val display = crdtmap.triangulation.getDisplay(poly);
  val b = poly.bounds;
  println(b);
  
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

/*

import scala.swing._
import java.awt.{ Graphics2D, Color }

object MainMap extends SimpleSwingApplication {
  val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps);
  val map = DiskMap.fromPoly(poly);
  
  println("end newton");
  println(map.C);
  println(map.A);
  for(id <- 0 until poly.n){
    print(map.pre_apply(id) * map.C + map.A);
    print(", ");
  }
    print("\n");
    
  for(id <- 0 until poly.n){
    print(map.angles(id));
    print(", ");
  }
    print("\n");
    
  val display = map.getDisplay;
  val b = poly.bounds;
  println(b);
  
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

/*
import spire.math.Complex
import spire.implicits._


class DiskMap(val angles : IndexedSeq[Double], val pre_angles : Array[Double], var C : Complex[Double], var A : Complex[Double]) {
	
	def derivative(z : Complex[Double]) :  Complex[Double] = pre_derivative(z) * C;
	def pre_derivative(z : Complex[Double]) :  Complex[Double] = {
	  var out : Complex[Double] = 1;
	  
	  for(id <- 0 until angles.length){
	    out *= (1.0 - z / Complex.polar(1.0, pre_angles(id)))**(angles(id)/Math.PI - 1.0);
	  }
	  return out;
	}
	
	def pre_derivative_without(z : Complex[Double], without : Int) : Complex[Double] = {
	  var out : Complex[Double] = 1;
	  
	  for(id <- 0 until angles.length){
	    if(id != without) {
	      val a = Complex.polar(1.0, pre_angles(id));
	      if(a  == z){
	        println("oups");
	        println(id);
	        println(pre_angles(id));
	        println(pre_angles(without));
	      }
	      out *= (1 - z / a)**(angles(id)/Math.PI - 1.0);
	    }
	  }
	  return out;
	}
	
	def pre_apply(a : Complex[Double], b : Complex[Double]) : Complex[Double] = {
	  if(a == b) return 0;
	  
	  var out : Complex[Double] = 0;
	  val d : Complex[Double] = b - a;
	  
	  val steps : Int = Math.max(1, (d.abs * 500).toInt);
	  val d_phi : Complex[Double] = d / steps;
	  
	  for(id <- 0 until steps){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  return out;
	}
	
	def apply(idA : Int, idB : Int) : Complex[Double] = pre_apply(idA, idB) * C + A;
	def pre_apply(idA : Int, idB : Int) : Complex[Double] = {
	  if(idA == idB) return 0;
	  
	  val a = Complex.polar(1.0, pre_angles(idA));
	  val b = Complex.polar(1.0, pre_angles(idB));
	  
	  var out : Complex[Double] = 0;
	  val d : Complex[Double] = b - a;
	  
	  val steps : Int = Math.max(2, (d.abs * 500).toInt);
	  val d_phi : Complex[Double] = d / steps;
	  
	  out += {
		  val z1 = a + d_phi;
		  val f0 = pre_derivative_without(a, idA);
		  val f1 = pre_derivative_without(z1, idA);
		  val r = angles(idA)/Math.PI;
		  -a * ((1 - z1 / a)**r) * f1 / r - a * a * (f1 - f0) * ((1 - z1 / a)**(r+1)) / (d_phi * r * (r+1));
	  }
	  // f'' = (1 - z / a)**(r - 1)
	  // f' = - a * (1 - z / a)**r / r
	  // f = a * a * (1 - z / a)**(r + 1) / (r* (r+1))
	  
	  for(id <- 1 until steps-1){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  
	  out += {
		  val z0 = b - d_phi;
		  val f0 = pre_derivative_without(z0, idB);
		  val f1 = pre_derivative_without(b, idB);
		  val r = angles(idB)/Math.PI;
		  b * ((1 - z0 / b)**r) * f0 / r + b * b * (f1 - f0) * ((1 - z0 / b)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  return out;
	}
	
	def apply(id : Int) : Complex[Double] = pre_apply(id) * C + A;
	def pre_apply(id : Int) : Complex[Double] = pre_apply(0.0, id);
	def apply(a : Complex[Double], idB : Int) : Complex[Double] = pre_apply(a, idB) * C;
	def apply(idA : Int, b : Complex[Double]) : Complex[Double] = pre_apply(idA, b) * C;
	def pre_apply(a : Complex[Double], idB : Int) : Complex[Double] = {
	  val b = Complex.polar(1.0, pre_angles(idB));
	  val d = b - a;
	  
	  var out : Complex[Double] = 0;
	  
	  val steps : Int = Math.max(2, (d.abs * 500).toInt);
	  val d_phi : Complex[Double] = d / steps.toDouble;
	  
	  for(id <- 0 until steps-1){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  out += {
		  val z0 = b - d_phi;
		  val f0 = pre_derivative_without(z0, idB);
		  val f1 = pre_derivative_without(b, idB);
		  val r = angles(idB)/Math.PI;
		  b * ((1 - z0 / b)**r) * f0 / r + b * b * (f1 - f0) * ((1 - z0 / b)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  return out;
	}
	
	def pre_apply(idA : Int, b : Complex[Double]) : Complex[Double] = {
	  val a = Complex.polar(1.0, pre_angles(idA));
	  val d = b - a;
	  
	  var out : Complex[Double] = 0;
	  
	  val steps : Int = Math.max(2, (d.abs * 500).toInt);
	  val d_phi : Complex[Double] = d / steps;
	  
	  out += {
		  val z1 = a + d_phi;
		  val f0 = pre_derivative_without(a, idA);
		  val f1 = pre_derivative_without(z1, idA);
		  val r = angles(idA)/Math.PI;
		  -a * ((1 - z1 / a)**r) * f1 / r - a * a * (f1 - f0) * ((1 - z1 / a)**(r+1)) / (d_phi * r * (r+1));
	  }
	  
	  for(id <- 1 until steps){
	    val phi0 = a + d * ((id.toDouble + 0.25) / steps.toDouble);
	    val phi1 = a + d * ((id.toDouble + 0.5) / steps.toDouble);
	    val phi2 = a + d * ((id.toDouble + 0.75) / steps.toDouble);
	    
	    out += (pre_derivative(phi0)*2 - pre_derivative(phi1) + pre_derivative(phi2)*2) * d_phi / 3;
	  }
	  
	  return out;
	}
	def apply(z : Complex[Double]) : Complex[Double] = pre_apply(z) * C + A;
	def pre_apply(z : Complex[Double]) : Complex[Double] = pre_apply(0.0, z);
	def map_apply(zs : Array[Complex[Double]]){
		// here we assume a path and integrate gradually
		// replacing with the result

		var last_w = apply(zs(0));
		var last_z = zs(0);
		zs(0) = last_w;
		
		for(id <- 1 until zs.length){
			val z = zs(id);
			
			val cuts = halfRule(last_z, z, 1);
			if(cuts.length == 0){
			    val z0 = (z + last_z * 3) / 4;
			    val z1 = (z + last_z) / 2;
			    val z2 = (z * 3 + last_z) / 4;
			    last_w += (pre_derivative(z0)*2 - pre_derivative(z1) + pre_derivative(z2)*2) * C * (z - last_z) / 3;
			    
			}else{
			  var lastID = cuts.head;
			  var last = Complex.polar(1.0, pre_angles(cuts.head));
			  
			  if(last != last_z) last_w += apply(last_z, lastID);
			    
			  for(nextID <- cuts.tail){
			    last_w += apply(lastID, nextID);
			    lastID = nextID;
			  }
			  
			  //lastID = cuts.last;
			  last = Complex.polar(1.0, pre_angles(lastID));
			  if(last != z) last_w += apply(lastID, z);
			  
			}
			zs(id) = last_w;
			last_z = z;
		}
	}
	
	def halfRule(a : Complex[Double], b : Complex[Double], steps : Int) : Seq[Int] = {
	  val border = 0.5 / steps.toDouble;
	  return pre_angles.view.map((t : Double) => (Complex.polar(1.0, t)-a) / (b - a))
	  			.zipWithIndex.collect({ case (Complex(x, y), id)
	  				if (-border <= x && x <= 1.0+border && -border <= y && y <= border)
	  				  => (x, id)
	  			}).sortBy(_._1).map(_._2);
	}
	
	// inversing the map
	def inverse(w : Complex[Double]) : Complex[Double] = inverse_newton(w, inverse_euler(w));
	def inverse_euler(w : Complex[Double]) : Complex[Double] = {
	   var last_z : Complex[Double] = 0;
	   var last_w = apply(last_z);
	   
	   val delta = w - last_w;
	   val steps : Int = 200;
	   val h : Complex[Double] = delta / steps;
	   
	   for(id <- 0 until steps){
	     last_z += h / derivative(last_z);
	   }
	   return last_z;
	}
	def inverse_newton(w : Complex[Double], guess : Complex[Double]) : Complex[Double] = {
	  var z = guess;
	  
	  val err = 0.001;
	  val sqr_err = err**2;
	  
	  for(id <- 0 until 300){
	    val d = apply(z) - w;
	    if(d.norm < sqr_err) return z;
	    
	    z -= d / derivative(z);
	  }
	  
	  return z;
	}
	def map_inverse(ws : Array[Complex[Double]]){
	  if(ws.length == 0) return;
	  
	  var last_z = inverse(ws(0));
	  ws(0) = last_z;
	  for(id <- 1 until ws.length){
	    last_z = inverse_newton(ws(id), last_z);
	    ws(id) = last_z;
	  }
	}
	
	
	def getDisplay() : DisplayArray = {
	  val num_radii = 100;
	  val num_peris = 50;
	  val displays = new Array[Display](num_radii + num_peris);
	  
	  
	  val steps_radius = 300;
	  for(id <- 0 until num_radii){
	    val dir = Complex.polar(1.0, id.toDouble * 2 * Math.PI / num_radii.toDouble);
	    val steps = (0 to steps_radius) map (dir * _.toDouble / steps_radius.toDouble) toArray;
	    map_apply(steps);
	    displays(id) = new Path(steps);
	  }
	  
	  val steps_peri = 800;
	  for(id <- 0 until num_peris){
	    val r = (id.toDouble + 1.0) / num_peris.toDouble;
	    val steps = (0 to steps_peri) map ((i) => Complex.polar(r, i.toDouble * 2 * Math.PI / steps_peri.toDouble)) toArray;
	    map_apply(steps);
	    displays(num_radii + id) = new ClosedPath(steps);
	  }
	  
	  return new DisplayArray(displays);
	}
}

import jvx.numeric.PnMatrix;

object DiskMap{
  
  abstract class Condition{
    def apply(map : DiskMap, v0 : Complex[Double], v0_abs : Double) : List[Double];
  }
  case class NonInf(val idA : Int, val idB : Int, val prop : Double) extends Condition {
    def apply(map : DiskMap, v0 : Complex[Double], v0_abs : Double) = List(map.pre_apply(idA, idB).abs / v0_abs - prop);
  }
  case class Inf(val idA : Int, val idB : Int, val prop : Complex[Double]) extends Condition {
    def apply(map : DiskMap, v0 : Complex[Double], v0_abs : Double) = { val tmp = map.pre_apply(idA, idB) / v0 - prop; List(tmp.real, tmp.imag)}
  }
  
  
  def fromPoly(poly : Polygon) : DiskMap = {
    val n = poly.n;
    
    val pre_angles = new Array[Double](n);
    for(id <- 0 until n-3) pre_angles(id) = Math.PI*(id.toDouble+1) / (n-2).toDouble;
    pre_angles(n-3) = Math.PI * 1.0;
    pre_angles(n-2) = Math.PI * 1.5;
    pre_angles(n-1) = Math.PI * 2.0;
    
    val map = new DiskMap(poly.inner_angles, pre_angles, 1, 0);
    
    val w0_abs = (poly.points(1).unsafe_get - poly.points(0).unsafe_get).abs;
    val conditions = new Array[Condition](n-3);
    for(id <- 0 until n-3){
      conditions(id) = NonInf(id + 1, id + 2, (poly.points(id+2).unsafe_get - poly.points(id+1).unsafe_get).abs / w0_abs);
    }
    
    newton(poly, conditions, map);
    
    return map;
  }
  
  def condition(conditions : IndexedSeq[Condition], map : DiskMap) : IndexedSeq[Double] = {
    val v0 = map.pre_apply(0, 1);
    conditions.flatMap(_.apply(map, v0, v0.abs));
  }
  
  def unsafe_condition_lagrangien(conditions : IndexedSeq[Condition], map : DiskMap, phi : Array[Double]) : Array[Array[Double]] = {
    val n = map.pre_angles.length;
    val lagrangien = Array.ofDim[Double](n-3, n-3);
    
    val new_phi = new Array[Double](n-3);
    for(id <- 0 until n-3) new_phi(id) = phi(id);
    
    //val theta = new Array[Double](n);
    //for(id <- 0 until n) theta(id) = map.pre_angles(id);
    
    val cond0 = condition(conditions, map);//.toArray[Double];
    
    // discrete difference
    val d = 0.5;
    
    for(id <- 0 until n-3){
    	new_phi(id) = phi(id) + d;
    	//println(map.pre_angles(0));
    	phi2theta(map.pre_angles, new_phi);
    	/*
    	for(a <- map.pre_angles){
    	  print(a);
    	  print(", ");
    	}
    	print("\n");
    	*/
    	new_phi(id) = phi(id);
    	
    	val cond1 = condition(conditions, map);//.toArray[Double];
    	
    	for((idd, (c0, c1)) <- (0 until n-3) zip (cond0 zip cond1)){
    	  //println(c1);
    	  lagrangien(idd)(id) = (c1 - c0) / d;
    	}
    }
    
    return lagrangien;
  }
  
  def newton(poly : Polygon, conditions : IndexedSeq[Condition], map : DiskMap){
    val m = map.pre_angles.length-3;
    
    val err = 0.00001;
		
    val phi = new Array[Double](m);
    theta2phi(map.pre_angles, phi);
    
    for(id_steps <- 0 until 200){
      val temp = condition(conditions, map);
      println(temp);
      val cond = temp.toArray;
      
      if(temp.map(Math.abs(_)).max < err) {
        println(temp.map(Math.abs(_)).max);
        println(id_steps);
        
        map.C = (poly.points(1).unsafe_get - poly.points(0).unsafe_get) / map.pre_apply(0, 1);
        map.A = poly.points(0).unsafe_get - map.pre_apply(0) * map.C
        
        return;
      }
      
      val lagrangien = unsafe_condition_lagrangien(conditions, map, phi);
      
      for(as <- lagrangien){
		for(a <- as){
		  print(a);
		  print(", ");
		}
		print("\n");
      }
		print("\n");
		
      val indx = new Array[Int](m);
      
      PnMatrix.ludcmp(lagrangien, m, indx);
      PnMatrix.lubksb(lagrangien, m, indx, cond);
      
      for(id <- 0 until m) phi(id) -= cond(id);
      
      phi2theta(map.pre_angles, phi);
      
		for(a <- phi){
		  print(a);
		  print(", ");
		}
		print("\n");
      
		for(a <- map.pre_angles){
		  print(a);
		  print(", ");
		}
		print("\n");
    }
    
    map.C = (poly.points(1).unsafe_get - poly.points(0).unsafe_get) / map.pre_apply(0, 1);
    map.A = poly.points(0).unsafe_get - map.pre_apply(0) * map.C
  }
 
	def theta2phi(theta : Array[Double], phi : Array[Double]) {
		val n = theta.length;
		//double phi [] = new double[n-3];
		phi(0) = Math.log(theta(0)/(theta(1) - theta(0)));
		for(id <- 1 until (n-3)){
			phi(id) = Math.log((theta(id) - theta(id-1))/(theta(id+1) - theta(id)));
		}
	}

	def phi2theta(theta : Array[Double], phi : Array[Double]) {
		val n = phi.length;
		//double theta [] = new double[n+3];
		
		var  prod : Double = 1.0;
		var sum : Double = 1.0;
		
		for(id <- 0 until n){
			prod *= Math.exp(-phi(id));
			theta(id) = sum;
			sum += prod;
		}
		
		for(id <- 0 until n){
			theta(id) *= Math.PI / sum;
		}

		theta(n) = Math.PI;
		theta(n+1) = Math.PI * 3.0 / 2.0;
		theta(n+2) = Math.PI * 2.0;
	}
}

import scala.swing._
import java.awt.{ Graphics2D, Color }

object MainMap extends SimpleSwingApplication {
  val ps = Array(new Complex(-1.0, 1.0), new Complex(-1.0, 0.0), new Complex(0.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, -1.0), new Complex(1.0, 0.0), new Complex(2.0, 0.0), new Complex(2.0, 1.0), new Complex(1.0, 1.0), 
		  		new Complex(0.0, 2), new Complex(0.0, 1.0));
  val poly = new Polygon(ps.map(new NonInfP(_)));
  val map = DiskMap.fromPoly(poly);
  
  println("end newton");
  println(map.C);
  println(map.A);
  for(id <- 0 until poly.n){
    print(map.pre_apply(id) * map.C + map.A);
    print(", ");
  }
    print("\n");
    
  for(id <- 0 until poly.n){
    print(map.angles(id));
    print(", ");
  }
    print("\n");
    
  val display = map.getDisplay;
  val b = poly.bounds;
  println(b);
  
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