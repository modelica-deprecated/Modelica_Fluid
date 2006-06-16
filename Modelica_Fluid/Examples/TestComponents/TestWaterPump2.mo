model TestWaterPump2 
  import PC = Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics;
  function pumpConsChar2 = 
      Modelica_Fluid.SubClasses.Flowmachines.PumpCharacteristics.quadraticPower
      (
    q_nom={0,0.001,0.0015}, W_nom={550,650,800});
  extends TestWaterPumpDefault(Pump1(M=0.1,redeclare function 
        powerCharacteristic = 
          pumpConsChar2));
  
  annotation (Diagram);
  
end TestWaterPump2;
