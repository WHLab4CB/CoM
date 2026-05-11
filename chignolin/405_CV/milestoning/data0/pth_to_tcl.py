import torch
import numpy as np

def pth_to_tcl(pth_path, output_tcl_path):
    state_dict = torch.load(pth_path, map_location=torch.device('cpu'))

    with open(output_tcl_path, 'w') as f:
        for param_name, param_tensor in state_dict.items():
            param_np = param_tensor.numpy()

            if len(param_np.shape) == 1:
                tcl_value = "{ " + " ".join(f"{x:.8f}" for x in param_np) + " }"
            elif len(param_np.shape) == 2:
                tcl_value = "{\n"
                for row in param_np:
                    tcl_value += "    {" + " ".join(f"{x:.8f}" for x in row) + "}\n"
                tcl_value += "}"
            else:
                raise ValueError("Unsupported tensor dimension")

            tcl_var_name = param_name.replace('.','_')
            f.write(f"set {tcl_var_name} {tcl_value}\n\n")

if __name__ == "__main__":
    pth_to_tcl("best_model.pth", "network_params.tcl")
