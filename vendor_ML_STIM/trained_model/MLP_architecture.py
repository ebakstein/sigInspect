import torch
    
class MLP_STIM(torch.nn.Module):
    def __init__(self, input_size=9, output_size=1):
        super(MLP_STIM, self).__init__()
        self.fc1 = torch.nn.Linear(input_size, 7)
        self.fc2 = torch.nn.Linear(7, 4)
        self.fc3 = torch.nn.Linear(4, 4)
        self.fc4 = torch.nn.Linear(4, 2)
        self.fc5 = torch.nn.Linear(2, output_size)

    def forward(self, x):
        x = torch.tanh(self.fc1(x))
        x = torch.tanh(self.fc2(x))
        x = torch.tanh(self.fc3(x))
        x = torch.tanh(self.fc4(x))
        x = self.fc5(x)
        return x