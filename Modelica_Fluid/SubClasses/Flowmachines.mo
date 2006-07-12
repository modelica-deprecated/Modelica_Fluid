package Flowmachines 
  package PumpCharacteristics 
    function linearFlow "Linear flow characteristic" 
      extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseFlow;
      input SI.VolumeFlowRate q_nom[2] 
        "Volume flow rate for two operating points (single pump)";
      input SI.Height head_nom[2] "Pump head for two operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      /* Linear system to determine the coefficients:
  head_nom[1]*g = c[1] + q_nom[1]*c[2];
  head_nom[2]*g = c[1] + q_nom[2]*c[2];
  */
      Real c[2] = Modelica.Math.Matrices.solve([ones(2),q_nom],head_nom*g) 
        "Coefficients of linear head curve";
    algorithm 
      // Flow equation: head * g = q*c[1] + c[2];
      head := 1/g * (c[1] + q_flow*c[2]);
    end linearFlow;

    function quadraticFlow "Quadratic flow characteristic" 
      extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseFlow;
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)";
      input SI.Height head_nom[3] "Pump head for three operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  head_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3), q_nom, q_nom2],head_nom*g) 
        "Coefficients of quadratic head curve";
    algorithm 
      // Flow equation: head * g = c[1] + q_flow*c[2] + q_flow^2*c[3];
      head := 1/g * (c[1] + q_flow*c[2] + q_flow^2*c[3]);
    end quadraticFlow;

    function polynomialFlow "Polynomial flow characteristic" 
      extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseFlow;
      input SI.VolumeFlowRate q_nom[:] 
        "Volume flow rate for N operating points (single pump)";
      input SI.Height head_nom[:] "Pump head for N operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      Integer N = size(q_nom,1) "Number of nominal operating points";
      Real q_nom_pow[N,N] = {{q_nom[j]^(i-1) for j in 1:N} for i in 1:N} 
        "Rows: different operating points; columns: increasing powers";
      /* Linear system to determine the coefficients (example N=3):
  head_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[N] = Modelica.Math.Matrices.solve(q_nom_pow,head_nom*g) 
        "Coefficients of polynomial head curve";
    algorithm 
      // Flow equation (example N=3): head * g = c[1] + q_flow*c[2] + q_flow^2*c[3];
      // Note: the implementation is numerically efficient only for low values of Na
      head := 1/g * sum(q_flow^(i-1)*c[i] for i in 1:N);
    end polynomialFlow;

    function constantEfficiency "Constant efficiency characteristic" 
       extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseEfficiency;
       input Real eta_nom "Nominal efficiency";
    algorithm 
      eta := eta_nom;
    end constantEfficiency;

    function quadraticPower "Quadratic power consumption characteristic" 
      extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.basePower;
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)";
      input SI.Power W_nom[3] "Power consumption for three operating points";
    protected 
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  W_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  W_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  W_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3),q_nom,q_nom2],W_nom) 
        "Coefficients of quadratic power consumption curve";
    algorithm 
      consumption := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticPower;
  end PumpCharacteristics;
end Flowmachines;
