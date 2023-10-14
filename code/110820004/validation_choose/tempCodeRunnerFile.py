selected_columns = np.sort(np.concatenate([[0], random_normal + 1, random_normal + 2, random_tumor + normal_total + 1, random_tumor + normal_total + 2]))

# print(selected_columns)

# data = data_origin_df.iloc[:, selected_columns.tolist()]
# print(data)