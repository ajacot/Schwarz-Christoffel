package map

import spire.math.Complex
import spire.implicits._
import scala.Array.canBuildFrom

class Splitting (val poly : Polygon, val triangulation : Triangulation) {
	
	val protect = new Array[Boolean](poly.n);
	val id_cuts = new Array[Int](poly.n);
	val cuts = new Array[Seq[Complex[Double]]](poly.n);
	for(id <- 0 until poly.n){
	  cuts(id) = Seq();
	}
  
	def num = id_cuts(poly.n - 1);
	
	def cut_ends() = 
	  for((id0, id, id1) <- Ring.getTriples(0 until poly.n)){
	    if(poly.inner_angles(id) < Math.PI / 4.0){
	      protect(id) = true;
	      val x = poly.points(id);
	      val d0 = poly.points(id0) - x;
	      val d1 = poly.points(id1) - x;
	      val dd0 = d0.abs;
	      val dd1 = d1.abs;
	      val d = Math.min(dd0, dd1);
	      
	    }
	  }
	
	def add_cut(id : Int, x : Complex[Double]){
	  id_cuts.zipWithIndex.find(_._1 <= id) match {
	  	case Some((from, part)) => if(from == id) poly.points(part) else cuts(part)(id - from - 1);
	  	case None => return null;
	  }
	}
	def points = cuts.view.zipWithIndex.flatMap({case (cuts, part) => poly.points(part) +: cuts});
	
	def _id_cuts = cuts.scanLeft(0)((sum, next) => sum + next.length);
	/*
	def points(id : Int) : Complex[Double] = 
	  id_cuts.zipWithIndex.find(_._1 <= id) match {
	  	case Some((from, part)) => if(from == id) poly.points(part) else cuts(part)(id - from - 1);
	  	case None => return null;
	  }*/
	
}