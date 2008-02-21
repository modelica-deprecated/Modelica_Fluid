within FluidSandbox;
package ConnectionSemantics 
  "Classes to implement yet unsupported connection semantics" 
  extends Icons.VariantLibrary;
  
  model ConnectionSemantics 
    "Allows to implement unsupported connection semantics" 
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.ConnectionSemantics;
    
    annotation (
      defaultComponentName="semantics",
      Coordsys(extent=[-100,-20; 100,20], scale=0.05),
      Icon(Text(extent=[-3,0; 3,-64], string=" ")),
      Diagram);
    
  end ConnectionSemantics;
end ConnectionSemantics;
