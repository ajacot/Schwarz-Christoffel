package map

import spire.math.Complex
import spire.implicits._
import scala.swing._
import scala.Array.fallbackCanBuildFrom
import spire.math.Complex.doubleToComplex
import spire.math.Complex.intToComplex

@SerialVersionUID(100L)
class DoubleMap(val to : CRDTMap, val domain : Polygon, // domain is the polygon of the domain it is calculated beforehand
				val inner_angles0 : IndexedSeq[Double], // the alternative set of inner angles
				val Cs0 : Array[Complex[Double]], val As0 : Array[Complex[Double]] // the constants for each diagonal are saved here
		) extends Serializable{
  
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
  
  def map_apply(zs : Array[Complex[Double]], diag : Int = -1) {
    println("map_apply");
    val pre_vertices = Array.ofDim[Complex[Double]](inner_angles0.length-3, inner_angles0.length);
    
    for(id <- 0 until inner_angles0.length-3){
      to.triangulation.get_pre_angles(id, to.cross_ratios, pre_vertices(id));
    }
    
      
    var i_diag = nearest_diag(zs(0));
    
    val guess = if(diag > 0) 
    				(pre_vertices(i_diag)(diag) + pre_vertices(i_diag)((diag + 1) % inner_angles0.length)) / 2.0 
    			else null;
      
    var from_map = new DiskMap(inner_angles0, pre_vertices(i_diag), Cs0(i_diag), As0(i_diag));
    var to_map = new DiskMap(to.inner_angles, pre_vertices(i_diag), to.Cs(i_diag), to.As(i_diag));
    var last_w = if(guess != null) from_map.inverse(zs(0), guess) else from_map.inverse(zs(0));
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
  
  def getDisplay(num_lines : Int = 25, step_lines : Int = 200) : Display = {
    val (h_l, v_l) = domain.getGridE(num_lines, step_lines);
    return new DisplayArray(
        (h_l ++ v_l).map((line) =>
          new DisplayArray(line.map((part) => {
            val part_array = part._1.toArray;
            map_apply(part_array);
            new Path(part_array)
          }).toArray[Path])
        ).toArray[Display]
    );
  }

  def nearest_diag(z : Complex[Double]) : Int = to.triangulation.nearest_diag(z, domain);
}

object DoubleMap{
  def fromCRDT(to : CRDTMap, inner_angles0 : IndexedSeq[Double]) : DoubleMap = fromCRDT(to, inner_angles0, 0, 1, 0.0, 1.0)
  def fromCRDT(to : CRDTMap, inner_angles0 : IndexedSeq[Double], id0 : Int, id1 : Int, p0 : Complex[Double], p1 : Complex[Double]) : DoubleMap = {
    // the domain polygon will be moved so that the vertices at id0 and id1 will be p1 and p2
    val n = inner_angles0.length;
    val points = new Array[Complex[Double]](n);
    
    val Cs0 = new Array[Complex[Double]](n-3);
    val As0 = new Array[Complex[Double]](n-3);
    
    to.triangulation.setConstants(0, to.cross_ratios, inner_angles0, Cs0, As0, points);
    
    
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

import scala.util.Random
import java.io._

object MainDoubleMap extends SimpleSwingApplication {
  /*val ps = Array(new Complex(0.0, 0.0), new Complex(-1.0, 0.0), new Complex(0.0, -1.0), new Complex(1.0, 0.0), 
		  		new Complex(1.0, 2.0), new Complex(2.0, 2.0), new Complex(1.0, 3.0), new Complex(0.0, 2.0)
		  	);
  val inner_angles0 = Array(Math.PI / 2, Math.PI / 2, Math.PI, Math.PI, 
		  		Math.PI / 2, Math.PI / 2, Math.PI, Math.PI
		  	);
  val (poly, inner_angles) = Polygon.cut(new Polygon(ps), inner_angles0);
  */
  //val (poly, inner_angles) = Polygon.labyrinth(11, 11, 1.0, new Random(159)); // 59, 159
  val (poly, inner_angles) = Polygon.labyrinth(4, 4, 1.0, new Random(159)); // 59, 159
  
  //val (poly, inner_angles) = Polygon.cut(Polygon.spiral(17), Polygon.ancre_spiral(17));
  
  val map = CRDTMap.fromPoly(poly);
  
  //val inner_angles = map.triangulation.equi_inner_angles;
  /*
  // this allows to save an approximation
  val oos = new ObjectOutputStream(new FileOutputStream("D:/Arthur/uni/schwarz-christoffel/saves/full_labi2.crdtmap"))
  oos.writeObject(mmap)
  oos.close
*/
  // (3) read the object back in
  //val ois = new ObjectInputStream(new FileInputStream("D:/Arthur/uni/schwarz-christoffel/saves/full_labi2.crdtmap"));
  //val map : CRDTMap = ois.readObject().asInstanceOf[DoubleMap].to;
  
  //val inner_angles = Polygon.ancre_labi_split;
  /*
  val inner_angles = new Array[Double](poly.n);
  for(id <- 0 until poly.n){
    inner_angles(id) = Math.PI;
  }
  inner_angles(0)= Math.PI / 2.0;
  inner_angles(poly.n - 1) = Math.PI / 2.0;
  inner_angles(poly.n / 2 - 18) = Math.PI / 2.0;
  inner_angles(poly.n / 2 - 17) = Math.PI / 2.0;
  */
  
  val mmap = DoubleMap.fromCRDT(map, inner_angles.toArray, 0, 1, 0, 1);
  
  val display = mmap.getDisplay();
  val b = poly.bounds;// union mmap.domain.bounds;
 
  val bord : Double = 50;
  
  def top = new MainFrame {
    title = "hello";
    
    preferredSize = new Dimension(700, 700);
    
    contents = new Panel(){
      override def paintComponent(g : Graphics2D){
    	  g.setColor(new Color(0, 0, 0, 100));
    	  
    	  val scaling = Rect.getScaleInverse(b, new Rect(bord, bord, size.getWidth() - bord*2, size.getHeight() - bord*2));
    	  display.display(g, scaling);
    	  poly.getDisplay.display(g, scaling);
      }
    }
  }
}
