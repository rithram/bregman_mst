namespace bbtree {

   template<typename TDataType, class TBregmanDiv>
   int PartitionData<TDataType, TBregmanDiv>(std::vector<TDataType>& data)
   {
     
     // choose two starting means at random
     TDataType& mean_1;
     TDataType& mean_2;
     
     std::vector<int> assignments(data.size(), 0);
     
     bool converged = false;
     
     while(!converged)
     {
       
       // cooy the old ones for comparison
       std::vector<int> old_assignments = assignments;
       
       TDataType& point_1 = data[mean_1];
       TDataType& point_2 = data[mean_2];
       
       int num_1 = 0;
       int num_2 = 0;
       
       converged = true;
       
       for (int i = 0; i < data.size(); i++)
       {
         
         double dist_1 = TBregmanDiv::Divergence(data[i], point_1);
         double dist_2 = TBregmanDiv::Divergence(data[i], point_2);
         
         if (dist_1 < dist_2) {
           assignments[i] = 1;
           num_1++;
         }
         else {
           assignments[i] = 2;
           num_2++;
         }
         
         // Someone was assigned to a new mean
         if (old_assignments[i] != assignments[i]) 
         {
           converged = false;
         }
         
       } // loop over points to compute assignments
       
       // If we haven't converged, recompute the means
       // This assumes that TDataType lives in a group (and that x = 0 sets it to the identity)
       // If this isn't the case, can I still proceed?
       // The Bregman divergence has an x - y term, so this must be the case 
       if (!converged) { 
         
         mean_1 = 0;
         mean_2 = 0;

         for (int i = 0; i < data.size(); i++) {
         
           if (assignments[i] == 1) {
             mean_1 += data[i];
           }
           else {
             mean_2 += data[i];
           }
         
         }
       
         // now, normalize the new centroids
         // This assumes that division is meaningful for the data type as well.
         mean_1 /= num_1;
         mean_2 /= num_2;
       
       }
     
     } // loop until convergence
     
     // Now, partition the data vector based on the assignment values
     
     int left = 0;
     int right = data.size() - 1;
     while (left < right) {
       
       while (assignements[left] == 1)
       {
         left++
       }
       
       while(assignments[right] == 2)
       {
         right--;
       }
       
       // swap left and right
       TDataType temp = data[left];
       data[left] = data[right];
       data[right] = temp;
       
       assignments[left] = 1;
       assignments[right] = 2;
       
     }
     
     return left + 1;
     
   } // PartitionData


 } // namespace