function Verify_rd_lessthan_one(par)
   #Input:
   #   par: vector with parameters
   # Output:
   #   deterministic R0: deterministic squart reproductive number

   Rd0 = Compute_deterministic_R0(par)

   cond = Rd0 - 0.7
   test = sign(cond) == -1.0
   return test
end
