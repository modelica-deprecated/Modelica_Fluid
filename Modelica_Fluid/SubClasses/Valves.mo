package Valves 
package ValveCharacteristics "Functions for valve characteristics" 
    
  function linear "Linear characteristic" 
    extends Modelica_Fluid.Interfaces.Valves.ValveCharacteristics.baseFun;
  algorithm 
    rc := pos;
  end linear;
    
  function one "Constant characteristic" 
    extends Modelica_Fluid.Interfaces.Valves.ValveCharacteristics.baseFun;
  algorithm 
    rc := 1;
  end one;
    
  function quadratic "Quadratic characteristic" 
    extends Modelica_Fluid.Interfaces.Valves.ValveCharacteristics.baseFun;
  algorithm 
    rc := pos*pos;
  end quadratic;
    
  function equalPercentage "Equal percentage characteristic" 
    extends Modelica_Fluid.Interfaces.Valves.ValveCharacteristics.baseFun;
    input Real rangeability = 20 "Rangeability";
    input Real delta = 0.01;
  algorithm 
    rc := if pos > delta then rangeability^(pos-1) else 
            pos/delta*rangeability^(delta-1);
    annotation (Documentation(info="<html>
This characteristic is such that the relative change of the flow coefficient is proportional to the change in the stem position:
<p> d(rc)/d(pos) = k d(pos).
<p> The constant k is expressed in terms of the rangeability, i.e. the ratio between the maximum and the minimum useful flow coefficient:
<p> rangeability = exp(k) = rc(1.0)/rc(0.0).
<p> The theoretical characteristic has a non-zero opening when pos = 0; the implemented characteristic is modified so that the valve closes linearly when pos &lt delta.
</html>"));
  end equalPercentage;
    
end ValveCharacteristics;
end Valves;
