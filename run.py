import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("\033[1m" + "=" * 60)
print("OPTIMIZATION FRAMEWORK — Full Pipeline")
print("=" * 60 + "\033[0m")

# SECTION 1: Dataset Generation
from src import dataset_generation
ds_outputs = dataset_generation.run()

# SECTION 2: Cross-Validation
from src import cross_validation
cv_outputs = cross_validation.run(ds_outputs)

# Combine outputs for subsequent sections
all_outputs = {**ds_outputs, **cv_outputs}

# SECTION 3: Gap Detection
from src import gap_detection
all_outputs = gap_detection.run(all_outputs)

print("\n" + "\033[1m" + "=" * 60)
print("\033[32mPipeline complete!\033[0m")
print("=" * 60 + "\033[0m")